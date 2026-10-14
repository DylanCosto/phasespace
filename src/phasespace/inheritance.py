"""Phase evidence from parents' genotypes (transmission).

For a heterozygous child call with alleles (a0, a1) in genotype order, the
allele that came from the father is determined when only one assignment fits
the parents' genotypes: a0 from the father and a1 from the mother, or the
reverse. Two calls whose paternal alleles have indexes f1 and f2 lie on the
paternal copy together, so their blocks have relative orientation f1 XOR f2.
Sites where both parents are heterozygous are uninformative; impossible
combinations are reported as Mendelian conflicts and never used.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import pysam

from .inputs import _BASES, _number, _resolve_contig
from .models import PhaseSpaceError


def _usable(record: object, call: object, min_gq: float, min_dp: int, block: bool, reference: bool) -> bool:
    """A parent's genotype is trusted when it passes filters and GQ. Depth is
    checked when present: merged trio files often keep GQ but drop DP for
    genotypes filled in from gVCF reference blocks."""
    filters = set(record.filter.keys()) - ({"PASS", ".", "RefCall"} if reference else {"PASS", "."})
    ft = call.get("FT")
    failed = [item for item in (ft if isinstance(ft, (tuple, list)) else [ft]) if item not in (None, "", ".", "PASS")]
    gq = _number(call.get("GQ"))
    depth = _number(call.get("MIN_DP")) if block else None
    depth = depth if depth is not None else _number(call.get("DP"))
    return not filters and not failed and gq is not None and gq >= min_gq and (depth is None or depth >= min_dp)


def _genotype(reader: pysam.VariantFile, sample: str, contig: str, variant: dict, min_gq: float, min_dp: int):
    """A parent's allele sequences at a child call, or None when unknown."""
    start, end = variant["pos"], variant["pos"] + len(variant["ref"])
    try:
        records = list(reader.fetch(contig, start, end))
    except (ValueError, OSError) as exc:
        raise PhaseSpaceError(f"The parents' VCF needs a tabix/CSI index: {exc}") from exc
    for record in records:
        call = record.samples[sample]
        gt = call.get("GT")
        if gt is None or len(gt) != 2 or None in gt:
            continue
        alts = record.alts or ()
        block = all(alt in ("<NON_REF>", "<*>", ".") for alt in alts) and not any(gt)
        if record.start == start and (record.ref or "").upper() == variant["ref"]:
            alleles = tuple((record.alleles[index] or "").upper() for index in gt)
            if all(_BASES.fullmatch(allele) for allele in alleles) and _usable(record, call, min_gq, min_dp, False, not any(gt)):
                return alleles
            return None
        if block and record.start <= start and record.stop >= end and _usable(record, call, min_gq, min_dp, True, True):
            return (variant["ref"], variant["ref"])
    return None


def _paternal_index(child: tuple[str, str], father, mother) -> tuple[int | None, str]:
    """Index of the child allele inherited from the father, with a status."""
    fits = []
    for paternal in (0, 1):
        from_father, from_mother = child[paternal], child[1 - paternal]
        if (father is None or from_father in father) and (mother is None or from_mother in mother):
            fits.append(paternal)
    if not fits:
        return None, "conflict"
    if father is None and mother is None:
        return None, "unknown"
    if len(fits) == 2:
        return None, "uninformative"
    return fits[0], "informative"


def inheritance_evidence(result: dict, parents_path: str | Path, father: str, mother: str,
                         min_gq: float = 20, min_dp: int = 5) -> dict:
    hets = sorted((v for v in result["variants"] if v["gt"][0] != v["gt"][1]), key=lambda v: v["pos"])
    evidence = {"source": Path(parents_path).name, "father": father, "mother": mother, "variants": {},
                "pairs": [], "within_blocks": {"agree": 0, "disagree": 0}, "conflicts": 0}
    try:
        reader = pysam.VariantFile(str(parents_path))
    except (ValueError, OSError) as exc:
        raise PhaseSpaceError(f"Cannot open the parents' VCF: {exc}") from exc
    with reader:
        samples = set(reader.header.samples)
        missing = [name for name in (father, mother) if name not in samples]
        if missing:
            raise PhaseSpaceError(f"Parent sample{'s' if len(missing) > 1 else ''} {', '.join(missing)} not in {Path(parents_path).name}.")
        contig = _resolve_contig(result["transcript"]["contig"], reader.header.contigs) if reader.header.contigs else result["transcript"]["contig"]
        if contig is None:
            raise PhaseSpaceError(f"Contig {result['transcript']['contig']!r} is absent from the parents' VCF.")
        paternal: dict[str, int] = {}
        for variant in hets:
            child = tuple(variant["alleles"])
            father_gt = _genotype(reader, father, contig, variant, min_gq, min_dp)
            mother_gt = _genotype(reader, mother, contig, variant, min_gq, min_dp)
            index, status = _paternal_index(child, father_gt, mother_gt)
            evidence["variants"][variant["id"]] = {"status": status, "paternal_allele": index,
                                                   "father": "/".join(father_gt) if father_gt else None,
                                                   "mother": "/".join(mother_gt) if mother_gt else None}
            if status == "conflict":
                evidence["conflicts"] += 1
            if index is not None:
                paternal[variant["id"]] = index
    blocks: dict[tuple[str, str], dict] = {}
    block_of = {v["id"]: v.get("block") for v in hets}
    number = lambda block: (int(block.split("-")[-1]) if block.split("-")[-1].isdigit() else 0, block)
    for first, second in combinations([v["id"] for v in hets if v["id"] in paternal], 2):
        parity = paternal[first] ^ paternal[second]
        block_a, block_b = block_of[first], block_of[second]
        if block_a is None or block_b is None:
            continue
        if block_a == block_b:
            evidence["within_blocks"]["agree" if parity == 0 else "disagree"] += 1
            continue
        key = tuple(sorted((block_a, block_b), key=number))
        entry = blocks.setdefault(key, {"block_a": key[0], "block_b": key[1], "same": 0, "opposite": 0, "variant_pairs": []})
        entry["same" if parity == 0 else "opposite"] += 1
        entry["variant_pairs"].append({"a": first, "b": second, "same": int(parity == 0), "opposite": parity})
    evidence["pairs"] = sorted(blocks.values(), key=lambda entry: -(entry["same"] + entry["opposite"]))
    return evidence


def add_inheritance_evidence(result: dict, parents_path: str | Path, father: str, mother: str,
                             min_gq: float = 20, min_dp: int = 5) -> dict:
    """Attach transmission evidence and the verdict among arrangements it allows.

    A relationship is supported when every informative pair of calls agrees;
    a single informative pair is enough because transmission is not sampled.
    """
    from .reads import apply_supported
    evidence = inheritance_evidence(result, parents_path, father, mother, min_gq, min_dp)
    result["inheritance_evidence"] = apply_supported(result, evidence, min_count=1, min_fraction=1.0)
    return result
