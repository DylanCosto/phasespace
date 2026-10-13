"""Read one diploid coding locus without treating missing records as reference.

CDS coordinates are zero-based and half-open. Complete, codon-aligned CDS
annotations are required; exon phase bases are retained in the joined sequence.
User files are opened read-only. The reference FASTA must already be indexed.
"""
from __future__ import annotations

from collections.abc import Iterable
import gzip
import math
from pathlib import Path
import re
from urllib.parse import unquote

import pysam

from .models import Evidence, Locus, PhaseSpaceError, Transcript, Variant


_BASES = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


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


def _read_annotation(path: Path, transcript_id: str) -> tuple[str, str, str, tuple[tuple[int, int], ...]]:
    features: dict[str, dict[str, str]] = {}
    selected: list[tuple[str, int, int, str, int, dict[str, str]]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("##FASTA"):
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise PhaseSpaceError(f"Annotation line {line_number} must contain nine tab-separated fields.")
                attrs = _attributes(fields[8])
                # Keep only metadata needed to resolve a transcript's gene;
                # exon/CDS feature IDs can otherwise dominate annotation memory.
                if "ID" in attrs and (fields[2] in ("gene", "pseudogene") or attrs["ID"] == transcript_id):
                    features[attrs["ID"]] = attrs
                if fields[2] in ("transcript", "mRNA") and attrs.get("transcript_id") == transcript_id:
                    features[transcript_id] = attrs
                parents = attrs.get("transcript_id", attrs.get("Parent", "")).split(",")
                if transcript_id in parents and fields[2].lower() == "selenocysteine":
                    raise PhaseSpaceError(f"Transcript {transcript_id!r} has a Selenocysteine feature; translation exceptions are unsupported under the standard genetic code model.")
                if fields[2] != "CDS" or transcript_id not in parents:
                    continue
                try:
                    start, end = int(fields[3]) - 1, int(fields[4])
                    phase = int(fields[7])
                except ValueError as exc:
                    raise PhaseSpaceError(f"CDS line {line_number} has invalid coordinates or missing numeric phase.") from exc
                if start < 0 or end <= start or phase not in (0, 1, 2):
                    raise PhaseSpaceError(f"CDS line {line_number} has invalid coordinates or phase.")
                if fields[6] not in ("+", "-"):
                    raise PhaseSpaceError(f"CDS line {line_number} must specify + or - strand.")
                selected.append((fields[0], start, end, fields[6], phase, attrs))
    except (OSError, UnicodeError) as exc:
        raise PhaseSpaceError(f"Cannot read annotation {path}: {exc}") from exc
    if not selected:
        raise PhaseSpaceError(f"No CDS features found for transcript {transcript_id!r}; use the exact transcript_id or GFF3 Parent identifier.")
    if len({row[0] for row in selected}) != 1 or len({row[3] for row in selected}) != 1:
        raise PhaseSpaceError("A transcript must have CDS features on one contig and one strand.")
    selected.sort(key=lambda row: row[1])
    if any(left[2] > right[1] for left, right in zip(selected, selected[1:])):
        raise PhaseSpaceError("Overlapping or duplicate CDS annotation intervals are unsupported.")
    strand = selected[0][3]
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
    attrs = ordered[0][5]
    transcript_attrs = features.get(transcript_id, {})
    for attributes in [transcript_attrs, *(row[5] for row in selected)]:
        translation_table = attributes.get("transl_table")
        if translation_table is not None and translation_table != "1":
            raise PhaseSpaceError(f"Transcript {transcript_id!r} declares unsupported translation table {translation_table!r}; only the standard genetic code (table 1) is supported.")
        if "transl_except" in attributes:
            raise PhaseSpaceError(f"Transcript {transcript_id!r} declares transl_except; translation exceptions are unsupported under the standard genetic code model.")
    gene_id = attrs.get("gene_id") or attrs.get("gene") or transcript_attrs.get("Parent", "").split(",")[0]
    gene_attrs = features.get(gene_id, {})
    gene = attrs.get("gene_name") or gene_attrs.get("Name") or gene_id or transcript_id
    return selected[0][0], strand, gene, tuple((row[1], row[2]) for row in selected)


def _clips(start: int, end: int, intervals: tuple[tuple[int, int], ...]) -> Iterable[tuple[int, int]]:
    for left, right in intervals:
        if start < right and end > left:
            yield max(start, left), min(end, right)


def _number(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        if len(value) != 1:
            return None
        value = value[0]
    try:
        result = float(value)
        return result if result >= 0 and result < float("inf") else None
    except (ValueError, TypeError):
        return None


def _quality(record: object, call: object, min_gq: float, min_dp: int, block: bool) -> tuple[list[str], float | None, int | None]:
    reasons = []
    filters = set(record.filter.keys())
    if filters - {"PASS", "."}:
        reasons.append("record FILTER=" + ";".join(sorted(filters)))
    ft = call.get("FT")
    ft_items = ft if isinstance(ft, (tuple, list)) else [ft]
    failed_ft = [str(item) for item in ft_items if item is not None and str(item) not in ("", ".", "PASS")]
    if failed_ft:
        reasons.append("sample FT=" + ";".join(failed_ft))
    gq = _number(call.get("GQ"))
    minimum_depth = _number(call.get("MIN_DP")) if block else None
    needs_minimum = block and record.stop > record.start + len(record.ref) and min_dp > 0
    if minimum_depth is not None:
        depth = minimum_depth
    elif needs_minimum:
        depth = None
    else:
        depth = _number(call.get("DP"))
    dp = int(depth) if depth is not None else None
    if gq is None:
        reasons.append("missing GQ")
    elif gq < min_gq:
        reasons.append(f"GQ below {min_gq:g}")
    if dp is None:
        reasons.append("missing MIN_DP for multi-base reference block" if needs_minimum else "missing DP/MIN_DP" if block else "missing DP")
    elif dp < min_dp:
        reasons.append(f"depth below {min_dp}")
    return reasons, gq, dp


def _project(start: int, end: int, ref: str, alleles: tuple[str, str], transcript: Transcript) -> tuple[int, int, str, tuple[str, str]]:
    offset = 0
    intervals = transcript.cds_intervals if transcript.strand == "+" else tuple(reversed(transcript.cds_intervals))
    for left, right in intervals:
        if left <= start and end <= right:
            for allele in alleles:
                if len(allele) <= len(ref):
                    continue
                prefix = 0
                while prefix < min(len(ref), len(allele)) and ref[prefix] == allele[prefix]:
                    prefix += 1
                suffix = 0
                while suffix < min(len(ref), len(allele)) - prefix and ref[len(ref) - suffix - 1] == allele[len(allele) - suffix - 1]:
                    suffix += 1
                if prefix + suffix == len(ref) and start + prefix in (left, right):
                    raise PhaseSpaceError(f"Insertion at {transcript.contig}:{start + 1} lies on a CDS boundary; its exon assignment is unsupported.")
            cds_start = offset + (start - left if transcript.strand == "+" else right - end)
            oriented_ref = ref if transcript.strand == "+" else _reverse_complement(ref)
            oriented_alleles = alleles if transcript.strand == "+" else tuple(_reverse_complement(allele) for allele in alleles)
            return cds_start, cds_start + len(ref), oriented_ref, oriented_alleles
        offset += right - left
    raise PhaseSpaceError(f"Variant at {transcript.contig}:{start + 1} crosses a CDS boundary. Only edits wholly inside one CDS segment are supported.")


def _deleted_tail(start: int, ref: str, alleles: Iterable[str]) -> tuple[int, int] | None:
    """Return the reference bases removed by the shortest called deletion allele."""
    tails = []
    for allele in alleles:
        if len(allele) >= len(ref):
            continue
        prefix = 0
        while prefix < len(allele) and ref[prefix] == allele[prefix]:
            prefix += 1
        tails.append(start + prefix)
    return (min(tails), start + len(ref)) if tails else None


def _without_deleted_tails(raw: list[Evidence], tails: list[tuple[int, int]]) -> list[Evidence]:
    """Drop reference calls over bases that a called deletion removes.

    gVCF writers (GATK HaplotypeCaller, DeepVariant) continue reference blocks
    through the bases deleted on one copy. The deletion record already carries
    that span as variant evidence, so these blocks are not conflicting calls.
    The anchor base is kept: a reference block covering it remains a conflict.
    """
    result = []
    for item in raw:
        pieces = [(item.start, item.end)]
        if item.status == "reference":
            for tail_start, tail_end in tails:
                pieces = [
                    piece
                    for left, right in pieces
                    for piece in ((left, min(right, tail_start)), (max(left, tail_end), right))
                    if piece[0] < piece[1]
                ]
        result.extend(Evidence(left, right, item.status, item.detail) for left, right in pieces)
    return result


def _merge_evidence(raw: list[Evidence], intervals: tuple[tuple[int, int], ...]) -> list[Evidence]:
    """Partition CDS, giving explicit uncertainty precedence over other records."""
    events: dict[int, list[tuple[str, str, int]]] = {}
    for left, right in intervals:
        events.setdefault(left, []).append(("cds", "", 1))
        events.setdefault(right, []).append(("cds", "", -1))
    for item in raw:
        events.setdefault(item.start, []).append((item.status, item.detail, 1))
        events.setdefault(item.end, []).append((item.status, item.detail, -1))
    boundaries = sorted(events)
    active: dict[str, dict[str, int]] = {"reference": {}, "variant": {}, "unknown": {}}
    cds_depth = 0
    result: list[Evidence] = []
    for start, end in zip(boundaries, boundaries[1:]):
        for status, detail, delta in events[start]:
            if status == "cds":
                cds_depth += delta
                continue
            bucket = active[status]
            bucket[detail] = bucket.get(detail, 0) + delta
            if bucket[detail] == 0:
                del bucket[detail]
        if not cds_depth:
            continue
        if active["variant"] and active["reference"]:
            raise PhaseSpaceError(f"Conflicting passing reference and variant genotypes overlap at genomic position {start + 1}. Resolve overlapping records before analysis.")
        status = next((name for name in ("unknown", "variant", "reference") if active[name]), "unknown")
        detail = "; ".join(sorted(item for item in active[status] if item)) if active[status] else "No usable genotype evidence in the input."
        if result and result[-1].end == start and result[-1].status == status and result[-1].detail == detail:
            previous = result.pop()
            result.append(Evidence(previous.start, end, status, detail))
        else:
            result.append(Evidence(start, end, status, detail))
    return result


def load_locus(vcf_path: str | Path, fasta_path: str | Path, annotation_path: str | Path,
               transcript_id: str, sample: str | None = None, min_gq: float = 20,
               min_dp: int = 5) -> Locus:
    """Load supported calls and explicit evidence for one annotated transcript.

    A missing quality field fails the quality threshold even if it is set to zero.
    Multi-base reference blocks require MIN_DP when min_dp is positive, since
    anchor DP cannot establish the minimum depth throughout the entire block.
    Phased records preserve their genotype order and PS. Phased calls without PS
    remain marked as such for the analysis layer's documented implicit block.
    """
    try:
        valid_thresholds = math.isfinite(min_gq) and min_gq >= 0 and math.isfinite(min_dp) and min_dp >= 0 and int(min_dp) == min_dp
    except (ValueError, TypeError, OverflowError):
        valid_thresholds = False
    if not valid_thresholds:
        raise PhaseSpaceError("Quality thresholds must be finite and nonnegative; minimum depth must be an integer.")
    contig, strand, gene, intervals = _read_annotation(Path(annotation_path), transcript_id)
    warnings: list[str] = []
    variants: list[Variant] = []
    raw_evidence: list[Evidence] = []
    called_spans: list[tuple[int, int]] = []
    deleted_tails: list[tuple[int, int]] = []
    if not Path(str(fasta_path) + ".fai").is_file():
        raise PhaseSpaceError("Reference FASTA requires an existing .fai index; index it with samtools faidx before running PhaseSpace.")
    try:
        fasta = pysam.FastaFile(str(fasta_path))
    except (ValueError, OSError) as exc:
        raise PhaseSpaceError(f"Cannot open indexed reference FASTA: {exc}") from exc
    with fasta:
        if contig not in fasta.references:
            raise PhaseSpaceError(f"Annotation contig {contig!r} is absent from the reference FASTA. Contig names must match exactly.")
        contig_length = fasta.get_reference_length(contig)
        if intervals[-1][1] > contig_length:
            raise PhaseSpaceError("CDS annotation extends beyond the reference contig.")
        reference_cds = "".join(fasta.fetch(contig, left, right).upper() for left, right in intervals)
        if not _BASES.fullmatch(reference_cds):
            raise PhaseSpaceError("The reference CDS contains unsupported bases; only A, C, G, T and N are supported.")
        if strand == "-":
            reference_cds = _reverse_complement(reference_cds)
        if "N" in reference_cds:
            warnings.append("Reference CDS contains N bases; affected codons have uncertain sequence consequences.")
        transcript = Transcript(transcript_id, gene, contig, strand, intervals, reference_cds)
        try:
            reader = pysam.VariantFile(str(vcf_path))
        except (ValueError, OSError) as exc:
            raise PhaseSpaceError(f"Cannot open VCF/gVCF/BCF: {exc}") from exc
        with reader:
            if reader.header.contigs and contig not in reader.header.contigs:
                raise PhaseSpaceError(f"Annotation contig {contig!r} is absent from the VCF header. Contig names must match exactly.")
            if contig in reader.header.contigs:
                declared_length = reader.header.contigs[contig].length
                if declared_length is not None and declared_length != contig_length:
                    raise PhaseSpaceError(f"Contig length mismatch for {contig!r}: VCF declares {declared_length}, reference FASTA has {contig_length}.")
            samples = list(reader.header.samples)
            if sample is None:
                if len(samples) != 1:
                    raise PhaseSpaceError("Choose a sample explicitly when the input does not contain exactly one sample.")
                sample = samples[0]
            if sample not in samples:
                raise PhaseSpaceError(f"Sample {sample!r} is absent from the input.")
            try:
                records = reader.fetch(contig, intervals[0][0], intervals[-1][1])
            except (ValueError, OSError):
                records = iter(reader)
                warnings.append("Input has no usable regional index; streamed records. BGZF VCF with a tabix/CSI index or indexed BCF enables regional loading.")
            try:
                for record in records:
                    if record.contig != contig:
                        continue
                    start = record.start
                    ref = (record.ref or "").upper()
                    ref_end = start + len(ref)
                    span_end = max(ref_end, record.stop)
                    clips = list(_clips(start, span_end, intervals))
                    if not clips:
                        continue
                    location = f"{contig}:{start + 1}"
                    if not _BASES.fullmatch(ref) or start < 0 or ref_end > contig_length:
                        raise PhaseSpaceError(f"Invalid REF allele or reference bounds at {location}.")
                    observed_ref = fasta.fetch(contig, start, ref_end).upper()
                    if ref != observed_ref:
                        raise PhaseSpaceError(f"REF mismatch at {location}: VCF has {ref}, reference FASTA has {observed_ref}. Check genome build and reference.")
                    call = record.samples[sample]
                    gt = call.get("GT")
                    called = [] if gt is None else [value for value in gt if value is not None]
                    if called and len(gt) != 2:
                        raise PhaseSpaceError(f"Non-diploid genotype at {location}; only diploid loci are supported.")
                    if any(value < 0 or value >= len(record.alleles or ()) for value in called):
                        raise PhaseSpaceError(f"Genotype references an absent allele at {location}.")
                    has_alt = any(value > 0 for value in called)
                    if has_alt:
                        called_sequences = [(record.alleles[value] or "").upper() for value in called]
                        if "*" in called_sequences:
                            # A spanning-deletion allele restates a deletion already
                            # called upstream. It adds no edit of its own.
                            if set(called_sequences) - {"*", ref}:
                                raise PhaseSpaceError(f"Unsupported called allele at {location}; a spanning deletion combined with another alternate allele cannot be resolved without read-level phasing.")
                            if not any(left <= start and ref_end <= right for left, right in deleted_tails):
                                raise PhaseSpaceError(f"Unsupported called allele at {location}; this spanning-deletion allele is not covered by a called upstream deletion.")
                            clips = list(_clips(start, ref_end, intervals))
                            reasons, _, _ = _quality(record, call, min_gq, min_dp, False)
                            if reasons:
                                detail = f"{location}: " + ", ".join(reasons)
                                warnings.append("Excluded call; " + detail + ".")
                                raw_evidence.extend(Evidence(left, right, "unknown", detail) for left, right in clips)
                            else:
                                raw_evidence.extend(Evidence(left, right, "variant", "Spanning deletion from an upstream call.") for left, right in clips)
                            continue
                        if any(not _BASES.fullmatch(allele) for allele in called_sequences):
                            raise PhaseSpaceError(f"Unsupported called allele at {location}; symbolic, spanning-deletion and breakend alleles cannot be translated.")
                        # END on a sequence variant does not extend its replacement span.
                        clips = list(_clips(start, ref_end, intervals))
                        if not clips:
                            continue
                        if any(start < end and ref_end > left for left, end in called_spans):
                            raise PhaseSpaceError(f"Overlapping called coding variants at {location} are unsupported; normalize and resolve the records first.")
                        called_spans.append((start, ref_end))
                        tail = _deleted_tail(start, ref, called_sequences)
                        if tail is not None:
                            deleted_tails.append(tail)
                        # Check unsupported boundary edits even when quality later fails.
                        known_alleles = tuple(called_sequences + [ref] * (2 - len(called_sequences)))
                        _project(start, ref_end, ref, known_alleles, transcript)
                    block_alts = record.alts or ()
                    block = not has_alt and all(alt in ("<NON_REF>", "<*>", ".") for alt in block_alts)
                    if not has_alt and not block:
                        # A 0/0 structural-variant record does not establish that
                        # every base in its END span is reference at other sites.
                        clips = list(_clips(start, ref_end, intervals))
                        if not clips:
                            continue
                    reasons, gq, dp = _quality(record, call, min_gq, min_dp, block)
                    if gt is None or len(gt) != 2 or any(value is None for value in gt):
                        reasons.insert(0, "missing or partial genotype")
                    if reasons:
                        detail = f"{location}: " + ", ".join(reasons)
                        warnings.append("Excluded call; " + detail + ".")
                        raw_evidence.extend(Evidence(left, right, "unknown", detail) for left, right in clips)
                        continue
                    if not has_alt:
                        raw_evidence.extend(Evidence(left, right, "reference", "Passing homozygous-reference genotype.") for left, right in clips)
                        continue
                    alleles = tuple(record.alleles[value].upper() for value in gt)
                    cds_start, cds_end, cds_ref, cds_alleles = _project(start, ref_end, ref, alleles, transcript)
                    if reference_cds[cds_start:cds_end] != cds_ref:
                        raise PhaseSpaceError(f"Internal CDS projection mismatch at {location}.")
                    ps = call.get("PS")
                    phase_set = str(ps) if call.phased and ps is not None else None
                    variant_id = f"{location}:{ref}>{','.join(alleles)}"
                    variants.append(Variant(variant_id, start, ref, alleles, tuple(gt), call.phased, phase_set,
                                            cds_start, cds_end, cds_ref, cds_alleles, gq, dp))
                    raw_evidence.extend(Evidence(left, right, "variant", "Passing variant genotype.") for left, right in clips)
            except (ValueError, OSError) as exc:
                if isinstance(exc, PhaseSpaceError):
                    raise
                raise PhaseSpaceError(f"Cannot interpret variant records: {exc}") from exc
    variants.sort(key=lambda variant: variant.cds_start)
    evidence = _merge_evidence(_without_deleted_tails(raw_evidence, deleted_tails), intervals)
    if any(item.status == "unknown" for item in evidence):
        warnings.append("Some CDS bases lack usable genotype evidence. An unchanged predicted copy is reference-identical only with respect to included calls; additional variants may be unobserved.")
    return Locus(transcript, sample, variants, evidence, list(dict.fromkeys(warnings)), Path(vcf_path).name)
