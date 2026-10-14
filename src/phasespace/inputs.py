"""Read one diploid coding locus without treating missing records as reference.

CDS coordinates are zero-based and half-open. Complete, codon-aligned CDS
annotations are required; exon phase bases are retained in the joined sequence.
User files are opened read-only. The reference FASTA must already be indexed.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
import gzip
import math
from pathlib import Path
import re
from urllib.parse import unquote

import pysam

from .annotation import TranscriptModel, _read_annotation, read_transcript  # noqa: F401 (re-exported)
from .models import Evidence, Locus, PhaseSpaceError, Transcript, Variant


_BASES = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


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


def _contig_aliases(name: str) -> list[str]:
    """The name itself, then the same contig with or without a ``chr`` prefix."""
    bare = name[3:] if name.lower().startswith("chr") else name
    names = [name, bare, "chr" + bare]
    if bare.upper() in ("M", "MT"):
        names += ["chrM", "MT", "chrMT", "M"]
    return list(dict.fromkeys(names))


def _resolve_contig(name: str, available: Iterable[str]) -> str | None:
    available = set(available)
    return next((alias for alias in _contig_aliases(name) if alias in available), None)


def check_thresholds(min_gq: float, min_dp: int) -> None:
    """Reject quality thresholds that cannot be applied; shared by the CLI and loader."""
    try:
        valid = math.isfinite(min_gq) and min_gq >= 0 and math.isfinite(min_dp) and min_dp >= 0 and int(min_dp) == min_dp
    except (ValueError, TypeError, OverflowError):
        valid = False
    if not valid:
        raise PhaseSpaceError("Quality thresholds must be finite and nonnegative; minimum depth must be an integer.")


def _quality(record: object, call: object, min_gq: float, min_dp: int, block: bool,
             reference_call: bool = False) -> tuple[list[str], float | None, int | None]:
    """Reasons a call cannot be used. DeepVariant marks confident homozygous-
    reference calls FILTER=RefCall; on a reference genotype that filter is not a
    failure, but GQ and depth thresholds still apply."""
    reasons = []
    filters = set(record.filter.keys())
    allowed = {"PASS", ".", "RefCall"} if reference_call else {"PASS", "."}
    if filters - allowed:
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


SPLICE_EXON_BASES, SPLICE_INTRON_BASES = 3, 8


def _splice_windows(intervals: tuple[tuple[int, int], ...]) -> list[tuple[int, int, int, str]]:
    """Splice-region windows at junctions between coding exons, as
    (start, end, junction, side): 3 exonic and 8 intronic bases per junction,
    following the Sequence Ontology splice_region_variant definition."""
    windows = []
    for (_, left_end), (right_start, _) in zip(intervals, intervals[1:]):
        windows += [(left_end - SPLICE_EXON_BASES, left_end, left_end, "exon"),
                    (left_end, left_end + SPLICE_INTRON_BASES, left_end, "intron"),
                    (right_start - SPLICE_INTRON_BASES, right_start, right_start, "intron"),
                    (right_start, right_start + SPLICE_EXON_BASES, right_start, "exon")]
    return windows


def _splice_hit(record: object, call: object, ref: str, windows: list[tuple[int, int, int, str]], location: str,
                min_gq: float, min_dp: int) -> dict | None:
    """Describe a passing alternate call that overlaps a splice-region window."""
    start, end = record.start, record.start + len(ref)
    hits = [(junction, side) for left, right, junction, side in windows if start < right and end > left]
    gt = call.get("GT")
    if not hits or gt is None or len(gt) != 2 or None in gt or not any(value > 0 for value in gt):
        return None
    alleles = [(record.alleles[value] or "").upper() for value in gt]
    if any(not _BASES.fullmatch(allele) for allele in alleles) or _quality(record, call, min_gq, min_dp, False)[0]:
        return None
    junction = min(hits, key=lambda hit: min(abs(start - hit[0]), abs(end - hit[0])))[0]
    # Bases from the junction: the last exonic and first intronic base are both 1.
    distance = junction - end + 1 if end <= junction else start - junction + 1 if start >= junction else 1
    return {"id": f"{location}:{ref}>{','.join(sorted(set(alleles) - {ref}))}", "pos": start, "ref": ref,
            "alleles": alleles, "gt": list(gt), "phased": bool(call.phased),
            "where": "exon and intron" if {side for _, side in hits} == {"exon", "intron"} else hits[0][1],
            "junction": junction, "distance": distance}


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

    Some gVCF writers continue reference blocks through the bases deleted on
    one copy. (DeepVariant 1.10 does not; it starts the next block after them.) The deletion record already carries
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


def _phased_calls(reader: pysam.VariantFile, path: Path, sample: str, contig: str,
                  region: tuple[int, int]) -> dict[tuple[int, str], tuple[tuple[str, str], object]]:
    """Phased heterozygous genotypes in a region, keyed by (start, REF)."""
    if reader.header.contigs:
        resolved = _resolve_contig(contig, reader.header.contigs)
        if resolved is None:
            raise PhaseSpaceError(f"Contig {contig!r} is absent from the phase file header.")
        contig = resolved
    try:
        records = reader.fetch(contig, *region)
    except (ValueError, OSError):
        records = iter(_open_variants(path, "phase file"))  # Unindexed: stream the whole file.
    phased: dict[tuple[int, str], tuple[tuple[str, str], object]] = {}
    try:
        for record in records:
            if record.contig != contig or not (record.start < region[1] and record.stop > region[0]):
                continue
            call = record.samples[sample]
            gt = call.get("GT")
            if not call.phased or gt is None or len(gt) != 2 or None in gt or gt[0] == gt[1]:
                continue
            alleles = tuple((record.alleles[index] or "").upper() for index in gt)
            phased[(record.start, (record.ref or "").upper())] = (alleles, call.get("PS"))
    except (ValueError, OSError) as exc:
        raise PhaseSpaceError(f"Cannot interpret phase file records: {exc}") from exc
    return phased


def _apply_phase(variants: list[Variant], phased: dict[tuple[int, str], tuple[tuple[str, str], object]],
                 label: str, warnings: list[str]) -> list[Variant]:
    """Take GT order and PS for heterozygous calls from a separate phased VCF.

    Phasing tools write variant-only VCFs; the main input keeps the gVCF's
    evidence. Calls match by position, REF and called allele sequences, so allele
    numbering may differ between the files. When a phase file is given it is the
    only source of phase: a heterozygous call it does not phase becomes
    unphased, so phase sets from the two files are never mixed.
    """
    result, applied, unmatched = [], 0, 0
    for variant in variants:
        if variant.gt[0] == variant.gt[1]:
            result.append(variant)
            continue
        match = phased.get((variant.pos, variant.ref))
        if match is None or sorted(match[0]) != sorted(variant.alleles):
            result.append(replace(variant, phased=False, phase_set=None))
            unmatched += 1
            continue
        order, phase_set = match
        if order != variant.alleles:
            variant = replace(variant, gt=variant.gt[::-1], alleles=variant.alleles[::-1], cds_alleles=variant.cds_alleles[::-1])
        result.append(replace(variant, phased=True, phase_set=None if phase_set is None else str(phase_set)))
        applied += 1
    total = applied + unmatched
    note = f"Phase taken from {label}: {applied} of {total} heterozygous coding call{'s' if total != 1 else ''} phased"
    warnings.append(note + (f"; {unmatched} not phased there {'is' if unmatched == 1 else 'are'} treated as unphased." if unmatched else "."))
    return result


def _open_variants(path: str | Path, what: str) -> pysam.VariantFile:
    try:
        return pysam.VariantFile(str(path))
    except (ValueError, OSError) as exc:
        raise PhaseSpaceError(f"Cannot open {what}: {exc}") from exc


class LocusReader:
    """Open the reference, calls and optional phase file once; load many transcripts.

    A missing quality field fails the quality threshold even if it is set to zero.
    Multi-base reference blocks require MIN_DP when min_dp is positive, since
    anchor DP cannot establish the minimum depth throughout the entire block.
    Phased records preserve their genotype order and PS. Phased calls without PS
    remain marked as such for the analysis layer's documented implicit block.
    With ``phase_from``, phase for heterozygous calls comes only from that VCF.
    Contig names may differ between inputs by a ``chr`` prefix (and chrM/MT);
    reported positions use the VCF's name.
    """

    def __init__(self, vcf_path: str | Path, fasta_path: str | Path, sample: str | None = None,
                 min_gq: float = 20, min_dp: int = 5, phase_from: str | Path | None = None):
        check_thresholds(min_gq, min_dp)
        self.vcf_path, self.min_gq, self.min_dp = Path(vcf_path), min_gq, min_dp
        if not Path(str(fasta_path) + ".fai").is_file():
            raise PhaseSpaceError("Reference FASTA requires an existing .fai index; index it with samtools faidx before running PhaseSpace.")
        try:
            self.fasta = pysam.FastaFile(str(fasta_path))
        except (ValueError, OSError) as exc:
            raise PhaseSpaceError(f"Cannot open indexed reference FASTA: {exc}") from exc
        self.vcf = _open_variants(vcf_path, "VCF/gVCF/BCF")
        samples = list(self.vcf.header.samples)
        if sample is None:
            if len(samples) != 1:
                self.close()
                raise PhaseSpaceError("Choose a sample explicitly when the input does not contain exactly one sample.")
            sample = samples[0]
        if sample not in samples:
            self.close()
            raise PhaseSpaceError(f"Sample {sample!r} is absent from the input.")
        self.sample = sample
        self.phase = self.phase_sample = None
        self.source = self.vcf_path.name
        if phase_from is not None:
            self.phase_path = Path(phase_from)
            self.phase = _open_variants(phase_from, "phase file")
            phase_samples = list(self.phase.header.samples)
            self.phase_sample = sample if sample in phase_samples else phase_samples[0] if len(phase_samples) == 1 else None
            if self.phase_sample is None:
                self.close()
                raise PhaseSpaceError(f"Sample {sample!r} is absent from the phase file, which has {len(phase_samples)} samples.")
            self.phase_label = f"{Path(phase_from).name} (sample {self.phase_sample})"
            self.source += f" + phase from {Path(phase_from).name}"

    def close(self) -> None:
        for handle in (getattr(self, "fasta", None), getattr(self, "vcf", None), getattr(self, "phase", None)):
            if handle is not None:
                handle.close()

    def __enter__(self) -> "LocusReader":
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    def _records(self, contig: str, start: int, end: int, warnings: list[str]):
        try:
            return self.vcf.fetch(contig, start, end)
        except (ValueError, OSError):
            warnings.append("Input has no usable regional index; streamed records. BGZF VCF with a tabix/CSI index or indexed BCF enables regional loading.")
            return iter(_open_variants(self.vcf_path, "VCF/gVCF/BCF"))

    def load(self, model: TranscriptModel) -> Locus:
        """Load supported calls and explicit evidence for one annotated transcript."""
        min_gq, min_dp, sample, fasta, reader = self.min_gq, self.min_dp, self.sample, self.fasta, self.vcf
        contig, strand, gene, intervals = model.contig, model.strand, model.gene, model.intervals
        warnings: list[str] = []
        variants: list[Variant] = []
        raw_evidence: list[Evidence] = []
        called_spans: list[tuple[int, int]] = []
        deleted_tails: list[tuple[int, int]] = []
        fasta_contig = _resolve_contig(contig, fasta.references)
        if fasta_contig is None:
            raise PhaseSpaceError(f"Annotation contig {contig!r} is absent from the reference FASTA (also tried {', '.join(repr(a) for a in _contig_aliases(contig)[1:])}).")
        contig_length = fasta.get_reference_length(fasta_contig)
        if intervals[-1][1] > contig_length:
            raise PhaseSpaceError("CDS annotation extends beyond the reference contig.")
        reference_cds = "".join(fasta.fetch(fasta_contig, left, right).upper() for left, right in intervals)
        if not _BASES.fullmatch(reference_cds):
            raise PhaseSpaceError("The reference CDS contains unsupported bases; only A, C, G, T and N are supported.")
        if strand == "-":
            reference_cds = _reverse_complement(reference_cds)
        if "N" in reference_cds:
            warnings.append("Reference CDS contains N bases; affected codons have uncertain sequence consequences.")
        vcf_contig = fasta_contig
        if reader.header.contigs:
            vcf_contig = _resolve_contig(contig, reader.header.contigs)
            if vcf_contig is None:
                raise PhaseSpaceError(f"Annotation contig {contig!r} is absent from the VCF header (also tried {', '.join(repr(a) for a in _contig_aliases(contig)[1:])}).")
            declared_length = reader.header.contigs[vcf_contig].length
            if declared_length is not None and declared_length != contig_length:
                raise PhaseSpaceError(f"Contig length mismatch for {vcf_contig!r}: VCF declares {declared_length}, reference FASTA has {contig_length}.")
        if len({contig, fasta_contig, vcf_contig}) > 1:
            warnings.append(f"Contig names differ between inputs (annotation {contig!r}, reference {fasta_contig!r}, VCF {vcf_contig!r}); they were matched by their chr prefix.")
        transcript = Transcript(model.id, gene, vcf_contig, strand, intervals, reference_cds)
        windows = _splice_windows(intervals)
        splice_region: list[dict] = []
        records = self._records(vcf_contig, intervals[0][0], intervals[-1][1], warnings)
        try:
            for record in records:
                if record.contig != vcf_contig:
                    continue
                start = record.start
                ref = (record.ref or "").upper()
                ref_end = start + len(ref)
                span_end = max(ref_end, record.stop)
                clips = list(_clips(start, span_end, intervals))
                if windows and any(start < right and start + len(ref) > left for left, right, _, _ in windows):
                    hit = _splice_hit(record, record.samples[sample], ref, windows, f"{vcf_contig}:{start + 1}", min_gq, min_dp)
                    if hit:
                        splice_region.append(hit)
                if not clips:
                    continue
                location = f"{vcf_contig}:{start + 1}"
                if not _BASES.fullmatch(ref) or start < 0 or ref_end > contig_length:
                    raise PhaseSpaceError(f"Invalid REF allele or reference bounds at {location}.")
                observed_ref = fasta.fetch(fasta_contig, start, ref_end).upper()
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
                reasons, gq, dp = _quality(record, call, min_gq, min_dp, block, reference_call=not has_alt)
                if gt is None or len(gt) != 2 or any(value is None for value in gt):
                    reasons.insert(0, "missing or partial genotype")
                if reasons:
                    detail = f"{location}: " + ", ".join(reasons)
                    warnings.append("Excluded call; " + detail + ".")
                    raw_evidence.extend(Evidence(left, right, "unknown", detail) for left, right in clips)
                    continue
                if not has_alt:
                    detail = "Homozygous-reference genotype (DeepVariant RefCall)." if "RefCall" in record.filter.keys() else "Passing homozygous-reference genotype."
                    raw_evidence.extend(Evidence(left, right, "reference", detail) for left, right in clips)
                    continue
                alleles = tuple(record.alleles[value].upper() for value in gt)
                cds_start, cds_end, cds_ref, cds_alleles = _project(start, ref_end, ref, alleles, transcript)
                if reference_cds[cds_start:cds_end] != cds_ref:
                    raise PhaseSpaceError(f"Internal CDS projection mismatch at {location}.")
                ps = call.get("PS")
                phase_set = str(ps) if call.phased and ps is not None else None
                # Name the call by its alternate alleles, independent of genotype order.
                variant_id = f"{location}:{ref}>{','.join(sorted(set(alleles) - {ref}))}"
                variants.append(Variant(variant_id, start, ref, alleles, tuple(gt), call.phased, phase_set,
                                        cds_start, cds_end, cds_ref, cds_alleles, gq, dp))
                raw_evidence.extend(Evidence(left, right, "variant", "Passing variant genotype.") for left, right in clips)
        except (ValueError, OSError) as exc:
            if isinstance(exc, PhaseSpaceError):
                raise
            raise PhaseSpaceError(f"Cannot interpret variant records: {exc}") from exc
        if self.phase is not None:
            phased = _phased_calls(self.phase, self.phase_path, self.phase_sample, vcf_contig, (intervals[0][0], intervals[-1][1]))
            variants = _apply_phase(variants, phased, self.phase_label, warnings)
        variants.sort(key=lambda variant: variant.cds_start)
        evidence = _merge_evidence(_without_deleted_tails(raw_evidence, deleted_tails), intervals)
        if any(item.status == "unknown" for item in evidence):
            warnings.append("Some CDS bases lack usable genotype evidence. An unchanged predicted copy is reference-identical only with respect to included calls; additional variants may be unobserved.")
        if splice_region:
            count = len(splice_region)
            warnings.append(f"{count} called variant{'s lie' if count != 1 else ' lies'} in a splice region (within {SPLICE_EXON_BASES} exonic or {SPLICE_INTRON_BASES} intronic bases of a coding-exon junction). Effects on splicing are not modelled.")
        return Locus(transcript, sample, variants, evidence, list(dict.fromkeys(warnings)), self.source, splice_region)


def load_locus(vcf_path: str | Path, fasta_path: str | Path, annotation_path: str | Path,
               transcript_id: str, sample: str | None = None, min_gq: float = 20,
               min_dp: int = 5, phase_from: str | Path | None = None) -> Locus:
    """Load supported calls and explicit evidence for one annotated transcript.

    See ``LocusReader`` for how calls, evidence, phase and contig names are handled.
    """
    check_thresholds(min_gq, min_dp)
    model = read_transcript(Path(annotation_path), transcript_id)
    with LocusReader(vcf_path, fasta_path, sample, min_gq, min_dp, phase_from) as reader:
        return reader.load(model)
