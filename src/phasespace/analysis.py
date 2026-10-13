"""Enumerate diploid phase hypotheses and their joint coding consequences.

This module deliberately does not infer phase probabilities or biological function.
Translation starts in the supplied annotated CDS frame, uses the standard genetic
code, and stops at the first stop codon. It does not model splicing, nonsense-
mediated decay, alternative initiation, mitochondrial translation, or expression.
"""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from itertools import combinations, product
import json

from .models import Locus, PhaseSpaceError, Variant


MAX_BLOCKS = 12
_DNA = frozenset("ACGTN")
# NCBI translation table 1; order is T,C,A,G for each codon position.
_BASES = "TCAG"
_AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
_GENETIC_CODE = {
    "".join(codon): amino
    for codon, amino in zip(product(_BASES, repeat=3), _AMINO_ACIDS)
}
# Every codon over the validated alphabet; any N-containing codon is unknown.
_CODONS = {
    "".join(codon): _GENETIC_CODE.get("".join(codon), "X")
    for codon in product(_BASES + "N", repeat=3)
}


def _translation(cds: str) -> tuple[str, int | None, str]:
    """Return expressed prefix, one-based first stop, and all complete codons."""
    codon = _CODONS.get
    full = "".join([codon(cds[i : i + 3], "X") for i in range(0, len(cds) - 2, 3)])
    stop = full.find("*")
    if stop < 0:
        return full, None, full
    return full[:stop], stop + 1, full


def _validate(locus: Locus, max_blocks: int) -> tuple[str, list[Variant]]:
    if isinstance(max_blocks, bool) or not isinstance(max_blocks, int):
        raise PhaseSpaceError("max_blocks must be an integer.")
    if not 0 <= max_blocks <= MAX_BLOCKS:
        raise PhaseSpaceError(f"max_blocks must be between 0 and {MAX_BLOCKS}.")
    reference = locus.transcript.reference_cds.upper()
    if not reference or len(reference) % 3:
        raise PhaseSpaceError("The reference CDS must be nonempty and divisible by three.")
    if set(reference) - _DNA:
        raise PhaseSpaceError("The reference CDS contains unsupported bases; use A, C, G, T or N.")
    if locus.transcript.strand not in ("+", "-"):
        raise PhaseSpaceError("The transcript strand must be '+' or '-'.")

    identifiers: set[str] = set()
    variants = sorted(locus.variants, key=lambda v: (v.cds_start, v.cds_end, v.id))
    for variant in variants:
        if not variant.id or variant.id in identifiers:
            raise PhaseSpaceError("Variant identifiers must be nonempty and unique.")
        identifiers.add(variant.id)
        if len(variant.gt) != 2 or len(variant.alleles) != 2 or len(variant.cds_alleles) != 2:
            raise PhaseSpaceError(f"{variant.id}: only complete diploid genotypes are supported.")
        if any(isinstance(g, bool) or not isinstance(g, int) or g < 0 for g in variant.gt):
            raise PhaseSpaceError(f"{variant.id}: genotype alleles must be nonnegative integers.")
        if not 0 <= variant.cds_start < variant.cds_end <= len(reference):
            raise PhaseSpaceError(f"{variant.id}: the coding interval is outside the reference CDS.")
        if variant.cds_ref.upper() != reference[variant.cds_start : variant.cds_end]:
            raise PhaseSpaceError(f"{variant.id}: the coding reference allele does not match the CDS.")
        for allele in (variant.ref, *variant.alleles, *variant.cds_alleles):
            if set(allele.upper()) - _DNA:
                raise PhaseSpaceError(f"{variant.id}: unsupported symbolic or ambiguous allele.")
        for index, gt in enumerate(variant.gt):
            if gt == 0 and (
                variant.alleles[index].upper() != variant.ref.upper()
                or variant.cds_alleles[index].upper() != variant.cds_ref.upper()
            ):
                raise PhaseSpaceError(f"{variant.id}: GT=0 must carry the reference allele.")
        if variant.gt[0] == variant.gt[1] and (
            variant.alleles[0].upper() != variant.alleles[1].upper()
            or variant.cds_alleles[0].upper() != variant.cds_alleles[1].upper()
        ):
            raise PhaseSpaceError(f"{variant.id}: identical genotype indexes have different sequences.")

    # Overlapping VCF records need normalization or a more general allele graph.
    # Reject the entire locus instead of dropping inconvenient phase hypotheses.
    previous: Variant | None = None
    for variant in variants:
        if all(a.upper() == variant.cds_ref.upper() for a in variant.cds_alleles):
            continue
        if previous is not None and variant.cds_start < previous.cds_end:
            raise PhaseSpaceError(
                f"Overlapping coding edits ({previous.id}, {variant.id}) are unsupported; "
                "normalize or combine these records before analysis."
            )
        previous = variant
    return reference, variants


def _blocks(variants: list[Variant]) -> tuple[list[dict], dict[str, str]]:
    blocks: list[dict] = []
    by_group: dict[tuple[str, str], dict] = {}
    membership: dict[str, str] = {}
    for variant in variants:
        if variant.gt[0] == variant.gt[1]:
            continue
        if not variant.phased:
            key = ("unphased", variant.id)
            label = f"Unphased {variant.id}"
        elif variant.phase_set is None:
            key = ("implicit", "")
            label = "Phased calls without PS (implicit phase set)"
        else:
            key = ("PS", str(variant.phase_set))
            label = f"Phase set {variant.phase_set}"
        if key not in by_group:
            block = {"id": f"block-{len(blocks) + 1}", "label": label, "variant_ids": []}
            blocks.append(block)
            by_group[key] = block
        block = by_group[key]
        block["variant_ids"].append(variant.id)
        membership[variant.id] = block["id"]
    return blocks, membership


def _edit_origins(start: int, ref: str, alt: str) -> list[int | None]:
    """Track homologous reference bases for stop-codon consequence comparisons."""
    if len(ref) == len(alt):
        return list(range(start, start + len(ref)))
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while (
        suffix < min(len(ref), len(alt)) - prefix
        and ref[len(ref) - suffix - 1] == alt[len(alt) - suffix - 1]
    ):
        suffix += 1
    return (
        list(range(start, start + prefix))
        + [None] * (len(alt) - prefix - suffix)
        + list(range(start + len(ref) - suffix, start + len(ref)))
    )


def _consequences(
    reference: str,
    cds: str,
    origins: list[int | None],
    deltas: list[int],
    reference_translation: tuple[str, int | None, str],
    translation: tuple[str, int | None, str],
) -> list[str]:
    if cds == reference:
        return ["reference"]
    ref_protein, ref_first_stop, ref_full = reference_translation
    protein, first_stop, full = translation
    consequences: list[str] = []

    if any(delta % 3 for delta in deltas):
        consequences.append("frameshift" if sum(deltas) % 3 else "frame_restored")
    if any(delta > 0 and delta % 3 == 0 for delta in deltas):
        consequences.append("inframe_insertion")
    if any(delta < 0 and delta % 3 == 0 for delta in deltas):
        consequences.append("inframe_deletion")

    # A retained terminal stop moves after an in-frame insertion/deletion. Its
    # amino-acid index alone therefore cannot identify a gained or lost stop.
    aligned: dict[int, tuple[int, str]]
    if not any(deltas):
        # Substitutions keep every base at its reference coordinate.
        aligned = {index: (index, amino) for index, amino in enumerate(full)}
    else:
        aligned = {}
        for index, amino in enumerate(full):
            origin = origins[index * 3]
            if (
                origin is not None
                and origin % 3 == 0
                and origins[index * 3 + 1] == origin + 1
                and origins[index * 3 + 2] == origin + 2
            ):
                aligned[origin // 3] = (index, amino)
    if first_stop is not None:
        bases = origins[(first_stop - 1) * 3 : first_stop * 3]
        origin = bases[0]
        retained_or_uncertain_stop = (
            origin is not None
            and origin % 3 == 0
            and bases == list(range(origin, origin + 3))
            and ref_full[origin // 3] in "*X"
        )
        if not retained_or_uncertain_stop:
            consequences.append("stop_gained")
    if ref_first_stop is not None and aligned.get(ref_first_stop - 1, (None, None))[1] not in ("*", "X"):
        consequences.append("stop_lost")
    if reference[:3] == "ATG" and cds[:3] != "ATG" and "N" not in cds[:3]:
        consequences.append("start_lost")

    if any(
        amino != ref_full[index] and amino not in "*X" and ref_full[index] not in "*X"
        and alt_index < len(protein) and index < len(ref_protein)
        for index, (alt_index, amino) in aligned.items()
    ):
        consequences.append("missense")
    if "X" in full or "X" in ref_full:
        consequences.append("uncertain_codon")
    if not any(deltas) and full == ref_full:
        consequences.append("synonymous")
    elif protein == ref_protein:
        consequences.append("protein_sequence_unchanged")
    # Complex indels can change unaligned codons without any simple substitution
    # label. Their in-frame/frame-restored labels and full translated sequence
    # are retained rather than pretending the unchanged net length is benign.
    return consequences


def _haplotype(
    reference: str,
    edits: list[tuple[Variant, str]],
    reference_translation: tuple[str, int | None, str],
) -> dict:
    cds = reference
    origins: list[int | None] = list(range(len(reference)))
    deltas: list[int] = []
    for variant, allele in sorted(edits, key=lambda edit: edit[0].cds_start, reverse=True):
        start, end = variant.cds_start, variant.cds_end
        cds = cds[:start] + allele + cds[end:]
        origins[start:end] = _edit_origins(start, variant.cds_ref.upper(), allele)
        deltas.append(len(allele) - (end - start))
    translation = _translation(cds)
    protein, first_stop, _ = translation
    return {
        "variant_ids": [variant.id for variant, _ in edits],
        "alleles": [allele for _, allele in edits],
        "cds": cds,
        "protein": protein,
        "consequences": _consequences(reference, cds, origins, deltas, reference_translation, translation),
        "first_stop": first_stop,
        "protein_length": len(protein),
    }


def _outcome_key(haplotypes: list[dict]) -> str:
    payload = json.dumps(sorted(h["cds"] for h in haplotypes), separators=(",", ":"))
    return "sha256:" + sha256(payload.encode("ascii")).hexdigest()


def _protein_outcome_key(haplotypes: list[dict]) -> str:
    """Group protein predictions, ignoring labels for sequence-preserving edits."""
    unchanged_labels = {"reference", "synonymous", "protein_sequence_unchanged"}
    copies = [
        json.dumps(
            [h["protein"], h["first_stop"], sorted(set(h["consequences"]) - unchanged_labels)],
            separators=(",", ":"),
        )
        for h in haplotypes
    ]
    payload = json.dumps(sorted(copies), separators=(",", ":"))
    return "sha256:" + sha256(payload.encode("ascii")).hexdigest()


def _informative_pairs(blocks: list[dict], scenarios: list[dict]) -> list[dict]:
    result = []
    for block_a, block_b in combinations(blocks, 2):
        a, b = block_a["id"], block_b["id"]
        outcomes: dict[str, set[int]] = {}
        for scenario in scenarios:
            parity = scenario["orientations"][a] ^ scenario["orientations"][b]
            outcomes.setdefault(scenario["protein_outcome_key"], set()).add(parity)
        # An outcome compatible with both parities cannot be ruled out by this
        # measurement. Count only outcome pairs with disjoint parity support.
        score = sum(values == {0} for values in outcomes.values()) * sum(
            values == {1} for values in outcomes.values()
        )
        result.append(
            {
                "block_a": a,
                "block_b": b,
                "score": score,
                "description": (
                    f"Same versus opposite block orientation separates {score} "
                    f"pair{'s' if score != 1 else ''} of distinct predicted protein outcomes; "
                    "this is not a probability or a biological severity score."
                ),
            }
        )
    return sorted(result, key=lambda pair: (-pair["score"], pair["block_a"], pair["block_b"]))


def analyze(locus: Locus, max_blocks: int = 8) -> dict:
    """Analyze every permitted phase arrangement; never silently truncate.

    Copy labels are arbitrary. Fixing the first block's orientation removes only
    global copy swaps; every relative arrangement remains represented. Phased
    calls without PS share one implicit block, as defined by the VCF convention.
    All coding outcomes are conditional on the supplied called variants and CDS.
    """
    reference, variants = _validate(locus, max_blocks)
    blocks, membership = _blocks(variants)
    if len(blocks) > max_blocks:
        raise PhaseSpaceError(
            f"This transcript has {len(blocks)} independent phase blocks, exceeding "
            f"max_blocks={max_blocks}. Select a smaller locus or increase the limit "
            f"up to {MAX_BLOCKS}; no scenarios were discarded."
        )

    reference_translation = _translation(reference)
    scenarios = []
    assignments = product((0, 1), repeat=max(0, len(blocks) - 1))
    for index, tail in enumerate(assignments, start=1):
        bits = (0, *tail) if blocks else ()
        orientations = {block["id"]: bit for block, bit in zip(blocks, bits)}
        copy_edits: list[list[tuple[Variant, str]]] = [[], []]
        for variant in variants:
            bit = orientations[membership[variant.id]] if variant.id in membership else 0
            for copy in (0, 1):
                allele = variant.cds_alleles[copy ^ bit].upper()
                if allele != variant.cds_ref.upper():
                    copy_edits[copy].append((variant, allele))
        haplotypes = [_haplotype(reference, edits, reference_translation) for edits in copy_edits]
        scenarios.append(
            {
                "id": f"scenario-{index}",
                "orientations": orientations,
                "haplotypes": haplotypes,
                "outcome_key": _outcome_key(haplotypes),
                "protein_outcome_key": _protein_outcome_key(haplotypes),
                "unaffected_copies": sum(h["cds"] == reference for h in haplotypes),
            }
        )

    warnings = list(locus.warnings)
    has_evidence_warning = any(
        "evidence" in warning.lower()
        and any(word in warning.lower() for word in ("lack", "unknown", "missing", "unobserved"))
        for warning in warnings
    )
    if any(e.status == "unknown" for e in locus.evidence) and not has_evidence_warning:
        warnings.append("Some CDS bases lack usable genotype evidence; additional variants may be unobserved.")
    if reference[-3:] not in ("TAA", "TAG", "TGA"):
        warnings.append("Reference CDS has no confirmed terminal stop; some annotations exclude the stop codon.")
    if any(v.phased and v.phase_set is None and v.gt[0] != v.gt[1] for v in variants):
        warnings.append(
            "Phased heterozygous calls without PS are treated as one implicit phase set. "
            "Confirm that the producing caller follows this VCF convention."
        )
    if any("N" in h["cds"] for s in scenarios for h in s["haplotypes"]):
        warnings.append(
            "N-containing codons translate to X; amino-acid consequences are uncertain "
            "and an unknown codon could contain an earlier stop."
        )
    ref_protein, ref_first_stop, _ = reference_translation
    return {
        "schema_version": 1,
        "sample": locus.sample,
        "transcript": asdict(locus.transcript),
        "variants": [dict(asdict(v), block=membership.get(v.id)) for v in variants],
        "evidence": [asdict(e) for e in locus.evidence],
        "warnings": list(dict.fromkeys(warnings)),
        "source": locus.source,
        "blocks": blocks,
        "scenarios": scenarios,
        "outcome_count": len({s["outcome_key"] for s in scenarios}),
        "protein_outcome_count": len({s["protein_outcome_key"] for s in scenarios}),
        "scenario_count": len(scenarios),
        "reference_protein": ref_protein,
        "reference_first_stop": ref_first_stop,
        "reference_protein_length": len(ref_protein),
        "informative_pairs": _informative_pairs(blocks, scenarios),
    }
