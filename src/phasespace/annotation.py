"""Find coding transcripts in GTF/GFF3 annotation.

Three entry points share one assembly step, so every path applies the same
checks (consistent CDS phase, joined GTF stop codons, standard genetic code):

- ``read_transcript`` finds one transcript by ID;
- ``select_transcript`` picks a gene's representative transcript;
- ``catalogue`` assembles one representative transcript per protein-coding gene
  in two passes over the file, for whole-genome scans.

Coordinates are zero-based and half-open.
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from urllib.parse import unquote

from .models import PhaseSpaceError

_VERSION = re.compile(r"\.\d+$")
_GTF_TAG = re.compile(r'\btag "([^"]+)"')
_GTF_TRANSCRIPT = re.compile(r'transcript_id "([^"]+)"')
_GFF_PARENT = re.compile(r"(?:^|;)Parent=([^;\n]+)")
# Representative transcript tags, best first. Matching ignores case, spaces and underscores.
_PREFERRED = (("maneselect", "MANE Select"), ("ensemblcanonical", "Ensembl canonical"), ("refseqselect", "RefSeq Select"))
_GENE_FEATURES = ("gene", "pseudogene")
_TRANSCRIPT_FEATURES = ("transcript", "mRNA")


@dataclass(frozen=True)
class TranscriptModel:
    id: str
    gene: str
    contig: str
    strand: str
    intervals: tuple[tuple[int, int], ...]
    selection: str = ""


def _attributes(value: str) -> dict[str, str]:
    if "=" in value.split(";", 1)[0]:
        return {
            key: unquote(item)
            for field in value.rstrip(";").split(";")
            if "=" in field
            for key, item in [field.strip().split("=", 1)]
        }
    result = {}
    for field in value.rstrip(";").split(";"):
        match = re.match(r'\s*(\S+)\s+"?([^";]*?)"?\s*$', field)
        if match:
            result[match.group(1)] = match.group(2)
    return result


def _tags(column: str, attrs: dict[str, str]) -> set[str]:
    """All transcript tags, normalised; GTF repeats ``tag``, GFF3 lists them."""
    tags = _GTF_TAG.findall(column) or attrs.get("tag", "").split(",")
    return {re.sub(r"[\s_]", "", tag).lower() for tag in tags if tag}


def _lines(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("##FASTA"):
                    return
                if line.strip() and not line.startswith("#"):
                    yield line_number, line
    except (OSError, UnicodeError) as exc:
        raise PhaseSpaceError(f"Cannot read annotation {path}: {exc}") from exc


def _fields(line: str, line_number: int) -> list[str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        raise PhaseSpaceError(f"Annotation line {line_number} must contain nine tab-separated fields.")
    return fields


def _transcript_names(attrs: dict[str, str]) -> set[str]:
    names = {attrs.get("ID"), attrs.get("transcript_id")}
    if attrs.get("transcript_id") and attrs.get("version"):
        names.add(f"{attrs['transcript_id']}.{attrs['version']}")
    names.discard(None)
    return names


def _assemble(transcript_id: str, rows: list[tuple[list[str], dict[str, str], int]],
              transcript_attrs: dict[str, str], gene_attrs_by_id: dict[str, dict[str, str]]) -> tuple[str, str, str, tuple[tuple[int, int], ...]]:
    """Build validated coding intervals from one transcript's CDS/stop_codon rows."""
    selected: list[tuple[str, int, int, str, int, dict[str, str]]] = []
    stops: list[tuple[str, int, int, str]] = []
    for fields, attrs, line_number in rows:
        if fields[2].lower() == "selenocysteine":
            raise PhaseSpaceError(f"Transcript {transcript_id!r} has a Selenocysteine feature; translation exceptions are unsupported under the standard genetic code model.")
        try:
            start, end = int(fields[3]) - 1, int(fields[4])
            phase = int(fields[7]) if fields[2] == "CDS" else 0
        except ValueError as exc:
            raise PhaseSpaceError(f"CDS line {line_number} has invalid coordinates or missing numeric phase.") from exc
        if start < 0 or end <= start or phase not in (0, 1, 2):
            raise PhaseSpaceError(f"CDS line {line_number} has invalid coordinates or phase.")
        if fields[6] not in ("+", "-"):
            raise PhaseSpaceError(f"CDS line {line_number} must specify + or - strand.")
        if fields[2] == "CDS":
            selected.append((fields[0], start, end, fields[6], phase, attrs))
        else:
            stops.append((fields[0], start, end, fields[6]))
    if not selected:
        raise PhaseSpaceError(f"No CDS features found for transcript {transcript_id!r}; use its transcript ID as written in the annotation (versioned or not), or its GFF3 ID.")
    if len({row[0] for row in selected}) != 1 or len({row[3] for row in selected}) != 1:
        raise PhaseSpaceError("A transcript must have CDS features on one contig and one strand.")
    selected.sort(key=lambda row: row[1])
    if any(left[2] > right[1] for left, right in zip(selected, selected[1:])):
        raise PhaseSpaceError("Overlapping or duplicate CDS annotation intervals are unsupported.")
    contig, strand = selected[0][0], selected[0][3]
    ordered = selected if strand == "+" else list(reversed(selected))
    total = 0
    for index, row in enumerate(ordered):
        expected = (3 - total % 3) % 3
        if row[4] != expected:
            reason = "Partial CDS models are unsupported: the first CDS phase must be 0." if index == 0 else f"Inconsistent CDS phase: expected {expected}, found {row[4]} at {row[0]}:{row[1] + 1}."
            raise PhaseSpaceError(reason)
        total += row[2] - row[1]
    if total % 3:
        raise PhaseSpaceError("The joined CDS length must be divisible by three; partial CDS models are unsupported.")
    intervals = [[row[1], row[2]] for row in selected]
    for stop_contig, left, right, stop_strand in sorted(stops, key=lambda stop: stop[1]):
        if stop_contig != contig or stop_strand != strand:
            raise PhaseSpaceError("A transcript's stop_codon must be on the same contig and strand as its CDS.")
        if any(start <= left and right <= end for start, end in intervals):
            continue  # Already inside the CDS, as in GFF3.
        adjacent = next((interval for interval in intervals if (interval[1] == left if strand == "+" else interval[0] == right)), None)
        if adjacent is not None:
            adjacent[0], adjacent[1] = min(adjacent[0], left), max(adjacent[1], right)
        else:
            intervals.append([left, right])  # A stop codon split across exons.
    intervals.sort()
    if any(left[1] > right[0] for left, right in zip(intervals, intervals[1:])):
        raise PhaseSpaceError("The annotated stop_codon overlaps the CDS in an unsupported way.")
    if sum(end - start for start, end in intervals) % 3:
        raise PhaseSpaceError("The joined CDS and stop codon length must be divisible by three.")
    attrs = ordered[0][5]
    for attributes in [transcript_attrs, *(row[5] for row in selected)]:
        translation_table = attributes.get("transl_table")
        if translation_table is not None and translation_table != "1":
            raise PhaseSpaceError(f"Transcript {transcript_id!r} declares unsupported translation table {translation_table!r}; only the standard genetic code (table 1) is supported.")
        if "transl_except" in attributes:
            raise PhaseSpaceError(f"Transcript {transcript_id!r} declares transl_except; translation exceptions are unsupported under the standard genetic code model.")
    gene_id = attrs.get("gene_id") or attrs.get("gene") or transcript_attrs.get("Parent", "").split(",")[0]
    gene_attrs = gene_attrs_by_id.get(gene_id, {})
    gene = attrs.get("gene_name") or transcript_attrs.get("gene_name") or gene_attrs.get("Name") or transcript_attrs.get("gene") or gene_id or transcript_id
    return contig, strand, gene, tuple((start, end) for start, end in intervals)


def read_transcript(path: Path, transcript_id: str) -> TranscriptModel:
    """Locate one transcript's coding intervals in a GTF or GFF3 file.

    The transcript may be named by its exact ID, by a GFF3 feature ID such as
    ``transcript:ENST…`` or ``rna-NM_…``, by ``transcript_id`` with or without
    a separate ``version`` attribute, or without a version when exactly one
    annotated version exists. GTF ``stop_codon`` features are joined to the CDS
    because GTF excludes the terminal stop codon from CDS features.
    """
    base = _VERSION.sub("", transcript_id)
    genes: dict[str, dict[str, str]] = {}
    transcript_attrs: dict[str, str] = {}
    aliases = {transcript_id}
    unversioned: dict[str, tuple[set[str], dict[str, str]]] = {}
    rows: list[tuple[list[str], dict[str, str], set[str], int]] = []
    for line_number, line in _lines(Path(path)):
        # Almost every line belongs to another transcript. Skip those before
        # splitting and parsing; gene lines are kept to resolve the gene name.
        if base not in line and "\tgene\t" not in line and "\tpseudogene\t" not in line:
            continue
        fields = _fields(line, line_number)
        attrs = _attributes(fields[8])
        if "ID" in attrs and fields[2] in _GENE_FEATURES:
            genes[attrs["ID"]] = attrs
        if fields[2] in _TRANSCRIPT_FEATURES:
            names = _transcript_names(attrs)
            if transcript_id in names:
                aliases.update(names)
                transcript_attrs = attrs
            elif base == transcript_id:
                for name in names:
                    if _VERSION.sub("", name) == transcript_id:
                        unversioned[name] = (names, attrs)
        if fields[2] in ("CDS", "stop_codon") or fields[2].lower() == "selenocysteine":
            parents = set(attrs.get("transcript_id", attrs.get("Parent", "")).split(","))
            if base in "".join(parents):
                rows.append((fields, attrs, parents, line_number))
    if not any(fields[2] == "CDS" and parents & aliases for fields, _, parents, _ in rows):
        # A versionless request matches a versioned ID only when unambiguous.
        versions = sorted({name for name in unversioned if not name.startswith(("transcript:", "rna-"))} or unversioned)
        if len(versions) > 1:
            raise PhaseSpaceError(f"Transcript {transcript_id!r} matches several annotated versions ({', '.join(versions)}); give the exact versioned ID.")
        if versions:
            names, transcript_attrs = unversioned[versions[0]]
            aliases.update(names)
    mine = [(fields, attrs, line_number) for fields, attrs, parents, line_number in rows if parents & aliases]
    contig, strand, gene, intervals = _assemble(transcript_id, mine, transcript_attrs, genes)
    return TranscriptModel(transcript_id, gene, contig, strand, intervals)


def _read_annotation(path: Path, transcript_id: str) -> tuple[str, str, str, tuple[tuple[int, int], ...]]:
    model = read_transcript(path, transcript_id)
    return model.contig, model.strand, model.gene, model.intervals


@dataclass
class _Candidate:
    id: str
    gene: str
    gene_id: str
    names: set[str]
    attrs: dict[str, str]
    tags: set[str]
    coding: bool
    contig: str = ""

    def rank(self) -> tuple[int, str] | None:
        for index, (tag, label) in enumerate(_PREFERRED):
            if tag in self.tags:
                return index, label
        return None


def _candidates(path: Path, wanted: set[str] | None) -> tuple[dict[str, dict[str, str]], list[_Candidate]]:
    """Transcript and gene features, optionally limited to lines naming a wanted gene."""
    genes: dict[str, dict[str, str]] = {}
    found: list[_Candidate] = []
    for line_number, line in _lines(Path(path)):
        if wanted is not None and not any(name in line for name in wanted):
            continue
        feature = line.split("\t", 3)[2] if line.count("\t") >= 8 else ""
        if feature not in _GENE_FEATURES and feature not in _TRANSCRIPT_FEATURES:
            continue
        fields = _fields(line, line_number)
        attrs = _attributes(fields[8])
        if feature in _GENE_FEATURES:
            if "ID" in attrs:
                genes[attrs["ID"]] = attrs
            continue
        biotype = attrs.get("transcript_type") or attrs.get("transcript_biotype") or attrs.get("biotype")
        names = _transcript_names(attrs)
        identifier = attrs.get("transcript_id") or attrs.get("ID")
        if identifier and attrs.get("version") and "." not in identifier:
            identifier = f"{identifier}.{attrs['version']}"
        parent = attrs.get("Parent", "").split(",")[0]
        found.append(_Candidate(identifier or "", attrs.get("gene_name") or attrs.get("gene") or "", attrs.get("gene_id") or parent,
                                names, attrs, _tags(fields[8], attrs), biotype in (None, "protein_coding", "mRNA"), fields[0]))
    for candidate in found:
        if not candidate.gene:
            candidate.gene = genes.get(candidate.gene_id, {}).get("Name", "")
    return genes, found


def _gene_matches(candidate: _Candidate, gene: str) -> bool:
    keys = {candidate.gene, candidate.gene_id, _VERSION.sub("", candidate.gene_id), candidate.gene_id.split(":")[-1], candidate.gene_id.removeprefix("gene-")}
    return gene in keys or _VERSION.sub("", gene) in keys


def select_transcript(path: Path, gene: str) -> tuple[str, str]:
    """Pick a gene's representative coding transcript: MANE Select, then
    Ensembl canonical, then RefSeq Select. Returns (transcript ID, reason)."""
    _, found = _candidates(Path(path), {gene, _VERSION.sub("", gene)})
    mine = [candidate for candidate in found if _gene_matches(candidate, gene) and candidate.coding]
    if not mine:
        raise PhaseSpaceError(f"No protein-coding transcript found for gene {gene!r}; check the gene symbol or ID used by this annotation.")
    ranked = sorted((candidate.rank(), candidate.id) for candidate in mine if candidate.rank() is not None)
    if not ranked:
        options = ", ".join(sorted(candidate.id for candidate in mine)[:10])
        raise PhaseSpaceError(f"Gene {gene!r} has no MANE Select, Ensembl canonical or RefSeq Select transcript in this annotation; choose one with --transcript ({options}).")
    (_, label), identifier = ranked[0]
    return identifier, f"{label} transcript {identifier} selected for gene {gene}."


def catalogue(path: Path, genes: set[str] | None = None) -> tuple[list[TranscriptModel], list[tuple[str, str, str]]]:
    """One representative coding transcript per gene, assembled in two passes.

    Returns the transcripts and a list of (gene, reason, contig) for genes
    skipped because no representative transcript exists or its model is
    unsupported. The contig is empty for requested genes that were not found.
    """
    gene_attrs, found = _candidates(Path(path), genes)
    best: dict[str, _Candidate] = {}
    skipped: list[tuple[str, str, str]] = []
    by_gene: dict[str, list[_Candidate]] = {}
    for candidate in found:
        if candidate.coding and candidate.id:
            key = candidate.gene or candidate.gene_id
            if genes is None or any(_gene_matches(candidate, gene) for gene in genes):
                by_gene.setdefault(key, []).append(candidate)
    for key, candidates in by_gene.items():
        ranked = sorted(((c.rank(), c.id, c) for c in candidates if c.rank() is not None), key=lambda item: (item[0], item[1]))
        if ranked:
            best[key] = ranked[0][2]
        else:
            skipped.append((key, "no MANE Select, Ensembl canonical or RefSeq Select transcript", candidates[0].contig))
    if genes is not None:
        matched = {gene for gene in genes if any(_gene_matches(c, gene) for c in best.values())}
        skipped.extend((gene, "gene not found", "") for gene in sorted(genes - matched) if not any(gene == entry[0] for entry in skipped))
    alias = {name: key for key, candidate in best.items() for name in candidate.names | {candidate.id}}
    rows: dict[str, list[tuple[list[str], dict[str, str], int]]] = {key: [] for key in best}
    for line_number, line in _lines(Path(path)):
        if "\tCDS\t" not in line and "\tstop_codon\t" not in line and "elenocysteine\t" not in line:
            continue
        match = _GTF_TRANSCRIPT.search(line)
        parents = [match.group(1)] if match else (_GFF_PARENT.search(line.rsplit("\t", 1)[-1]) or [None, ""])[1].split(",")
        owners = {alias[parent] for parent in parents if parent in alias}
        if owners:
            fields = _fields(line, line_number)
            attrs = _attributes(fields[8])
            for owner in owners:
                rows[owner].append((fields, attrs, line_number))
    models = []
    for key, candidate in best.items():
        try:
            contig, strand, gene, intervals = _assemble(candidate.id, rows[key], candidate.attrs, gene_attrs)
        except PhaseSpaceError as exc:
            skipped.append((key, str(exc), candidate.contig))
            continue
        models.append(TranscriptModel(candidate.id, candidate.gene or gene, contig, strand, intervals, candidate.rank()[1]))
    models.sort(key=lambda model: (model.contig, model.intervals[0][0]))
    return models, sorted(skipped)
