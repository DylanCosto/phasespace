"""Count reads that show which alleles of heterozygous calls travel together.

For two heterozygous calls whose alleles are numbered in genotype order, a read
(or read pair) carrying allele i1 at the first and i2 at the second came from
one gene copy, so their blocks have relative orientation i1 XOR i2 (0: same,
1: opposite). These are observed counts, not probabilities, and they never
change the input calls. Within a block the same test checks the input phase.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

import pysam

from .inputs import _resolve_contig
from .models import PhaseSpaceError

FLANK = 10
MIN_MAPQ = 20
MIN_BASEQ = 13


def _expected(fasta: pysam.FastaFile, contig: str, variant: dict, flank: int) -> tuple[int, int, list[str]]:
    start, end = variant["pos"], variant["pos"] + len(variant["ref"])
    left = fasta.fetch(contig, max(0, start - flank), start).upper()
    right = fasta.fetch(contig, end, end + flank).upper()
    return start - len(left), end + len(right), [left + allele + right for allele in variant["alleles"]]


def _type_read(read: pysam.AlignedSegment, variant: dict, positions: dict[int, int], windows: list) -> int | None:
    """Index (0 or 1, genotype order) of the allele this read shows, if clear."""
    sequence = read.query_sequence
    if sequence is None:
        return None
    if len(variant["ref"]) == 1 and all(len(allele) == 1 for allele in variant["alleles"]):
        query = positions.get(variant["pos"])
        if query is None or (read.query_qualities is not None and read.query_qualities[query] < MIN_BASEQ):
            return None
        base = sequence[query].upper()
        return variant["alleles"].index(base) if base in variant["alleles"] else None
    for window_start, window_end, expected in windows:
        left, right = positions.get(window_start - 1), positions.get(window_end)
        if left is None or right is None:
            continue
        observed = sequence[left + 1:right].upper()
        matches = [index for index, haplotype in enumerate(expected) if haplotype == observed]
        if len(matches) == 1:
            return matches[0]
    return None


def read_evidence(result: dict, read_paths: list[str | Path], fasta_path: str | Path,
                  min_mapq: int = MIN_MAPQ) -> dict:
    """Read support for relative orientations of the analysis's phase blocks."""
    contig = result["transcript"]["contig"]
    hets = [v for v in result["variants"] if v["gt"][0] != v["gt"][1] and v["alleles"][0] != v["alleles"][1]]
    hets.sort(key=lambda v: v["pos"])
    evidence = {"files": [Path(path).name for path in read_paths], "fragments": 0, "variants": {}, "pairs": [], "within_blocks": {"agree": 0, "disagree": 0}}
    if len(hets) < 1 or not read_paths:
        return evidence
    starts = [v["pos"] for v in hets]
    fragments: dict[str, dict[str, int | None]] = {}
    with pysam.FastaFile(str(fasta_path)) as fasta:
        fasta_contig = _resolve_contig(contig, fasta.references) or contig
        windows = {v["id"]: [_expected(fasta, fasta_contig, v, flank) for flank in (FLANK, 3)] for v in hets}
        for path in read_paths:
            try:
                reads = pysam.AlignmentFile(str(path), reference_filename=str(fasta_path))
            except (ValueError, OSError) as exc:
                raise PhaseSpaceError(f"Cannot open reads {path}: {exc}") from exc
            with reads:
                read_contig = _resolve_contig(contig, reads.references)
                if read_contig is None:
                    raise PhaseSpaceError(f"Contig {contig!r} is absent from {Path(path).name}.")
                try:
                    alignments = reads.fetch(read_contig, starts[0], hets[-1]["pos"] + len(hets[-1]["ref"]))
                except ValueError as exc:
                    raise PhaseSpaceError(f"Reads need a .bai/.crai/.csi index for regional access: {exc}") from exc
                for read in alignments:
                    if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail or read.mapping_quality < min_mapq:
                        continue
                    covered = hets[bisect_left(starts, read.reference_start - 50):bisect_right(starts, read.reference_end)]
                    covered = [v for v in covered if v["pos"] >= read.reference_start and v["pos"] + len(v["ref"]) <= read.reference_end]
                    if not covered:
                        continue
                    positions = {reference: query for query, reference in read.get_aligned_pairs(matches_only=True)}
                    calls = fragments.setdefault(f"{Path(path).name}:{read.query_name}", {})
                    for variant in covered:
                        allele = _type_read(read, variant, positions, windows[variant["id"]])
                        if allele is None:
                            continue
                        if variant["id"] in calls and calls[variant["id"]] != allele:
                            calls[variant["id"]] = None  # Mates or pieces disagree: drop this site.
                        else:
                            calls[variant["id"]] = allele
    per_variant = {v["id"]: [0, 0] for v in hets}
    pair_counts: dict[tuple[str, str], list[int]] = {}
    block_of = {v["id"]: v.get("block") for v in hets}
    order = {v["id"]: index for index, v in enumerate(hets)}
    blocks: dict[tuple[str, str], dict] = {}

    def block_key(first: str, second: str) -> tuple[str, str]:
        number = lambda block: (int(block.split("-")[-1]) if block.split("-")[-1].isdigit() else 0, block)
        return tuple(sorted((first, second), key=number))

    for calls in fragments.values():
        typed = sorted(((vid, allele) for vid, allele in calls.items() if allele is not None), key=lambda item: order[item[0]])
        if typed:
            evidence["fragments"] += 1
        votes: dict[tuple[str, str], list[int]] = {}
        for vid, allele in typed:
            per_variant[vid][allele] += 1
        for (first, a), (second, b) in combinations(typed, 2):
            pair_counts.setdefault((first, second), [0, 0])[a ^ b] += 1
            block_a, block_b = block_of[first], block_of[second]
            if block_a is None or block_b is None:
                continue
            if block_a == block_b:
                evidence["within_blocks"]["agree" if a == b else "disagree"] += 1
            else:
                votes.setdefault(block_key(block_a, block_b), [0, 0])[a ^ b] += 1
        # Each fragment counts once per block pair, by the majority of its variant pairs.
        for key, (same, opposite) in votes.items():
            entry = blocks.setdefault(key, {"block_a": key[0], "block_b": key[1], "same": 0, "opposite": 0, "variant_pairs": []})
            if same != opposite:
                entry["same" if same > opposite else "opposite"] += 1
    evidence["variants"] = {vid: {"allele_reads": counts} for vid, counts in per_variant.items()}
    for (first, second), (same, opposite) in pair_counts.items():
        block_a, block_b = block_of[first], block_of[second]
        if block_a is not None and block_b is not None and block_a != block_b:
            blocks[block_key(block_a, block_b)]["variant_pairs"].append({"a": first, "b": second, "same": same, "opposite": opposite})
    evidence["pairs"] = sorted(blocks.values(), key=lambda entry: -(entry["same"] + entry["opposite"]))
    return evidence


def supported_constraints(evidence: dict, min_reads: int = 3, min_fraction: float = 0.9) -> list[dict]:
    """Block relationships the reads clearly support: at least ``min_reads``
    informative fragments with at least ``min_fraction`` agreeing."""
    constraints = []
    for pair in evidence.get("pairs", []):
        total = pair["same"] + pair["opposite"]
        majority = max(pair["same"], pair["opposite"])
        if total >= min_reads and majority / total >= min_fraction:
            constraints.append({"a": pair["block_a"], "b": pair["block_b"], "parity": 0 if pair["same"] >= pair["opposite"] else 1})
    return constraints


def add_read_evidence(result: dict, read_paths: list[str | Path], fasta_path: str | Path,
                      min_reads: int = 3, min_fraction: float = 0.9) -> dict:
    """Attach read evidence, the relationships it clearly supports, and the
    verdict among arrangements consistent with them. Input calls are unchanged."""
    from .analysis import summarize
    evidence = read_evidence(result, read_paths, fasta_path)
    supported = supported_constraints(evidence, min_reads, min_fraction)
    compatible = [s for s in result["scenarios"]
                  if all(s["orientations"].get(c["a"], 0) ^ s["orientations"].get(c["b"], 0) == c["parity"] for c in supported)]
    evidence.update(supported=supported, thresholds={"min_reads": min_reads, "min_fraction": min_fraction},
                    compatible_arrangements=len(compatible),
                    verdict_with_reads=summarize({"scenarios": compatible})["verdict"] if compatible else "conflicting")
    result["read_evidence"] = evidence
    return result
