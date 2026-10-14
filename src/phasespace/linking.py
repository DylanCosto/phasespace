"""Link phase blocks by clear read or inheritance evidence before enumeration.

Opt-in (``--link-evidence``). Without it, evidence only filters arrangements
after enumeration, which is impossible for loci over the block cap. With it,
each clearly supported relationship (the same thresholds the report uses)
merges two blocks, strongest evidence first; a relationship that contradicts
stronger ones is skipped and reported. Calls in a merged block are rewritten as
one phase set, so the resulting arrangements assume that evidence.
"""
from __future__ import annotations

from dataclasses import asdict, replace

from .analysis import MAX_BLOCKS, _blocks, _validate
from .models import Locus


def _preview(locus: Locus) -> dict:
    """The parts of an analysis that evidence needs, without enumerating."""
    _, variants = _validate(locus, MAX_BLOCKS)
    blocks, membership = _blocks(variants)
    return {"transcript": asdict(locus.transcript), "blocks": blocks,
            "variants": [dict(asdict(v), block=membership.get(v.id)) for v in variants]}


def link_blocks(locus: Locus, *, reads: list | None = None, fasta: str | None = None,
                parents: tuple[str, str] | None = None, parents_from: str | None = None,
                min_gq: float = 20, min_dp: int = 5) -> tuple[Locus, list[dict]]:
    preview = _preview(locus)
    if len(preview["blocks"]) < 2:
        return locus, []
    edges = []  # (strength, block_a, block_b, parity, source)
    if reads:
        from .reads import read_evidence, supported_constraints
        evidence = read_evidence(preview, reads, fasta)
        totals = {(p["block_a"], p["block_b"]): p["same"] + p["opposite"] for p in evidence["pairs"]}
        edges += [(totals[(c["a"], c["b"])], c["a"], c["b"], c["parity"], "reads") for c in supported_constraints(evidence, 3, 0.9)]
    if parents:
        from .inheritance import inheritance_evidence
        from .reads import supported_constraints
        evidence = inheritance_evidence(preview, parents_from, *parents, min_gq, min_dp)
        totals = {(p["block_a"], p["block_b"]): p["same"] + p["opposite"] for p in evidence["pairs"]}
        edges += [(totals[(c["a"], c["b"])], c["a"], c["b"], c["parity"], "parents") for c in supported_constraints(evidence, 1, 1.0)]
    parent = {block["id"]: block["id"] for block in preview["blocks"]}
    flip = {block["id"]: 0 for block in preview["blocks"]}
    sources: dict[str, set[str]] = {}
    skipped = 0

    def find(block: str) -> tuple[str, int]:
        parity = 0
        while parent[block] != block:
            parity ^= flip[block]
            block = parent[block]
        return block, parity

    for _, a, b, parity, source in sorted(edges, key=lambda edge: -edge[0]):
        root_a, flip_a = find(a)
        root_b, flip_b = find(b)
        if root_a == root_b:
            skipped += flip_a ^ flip_b != parity
            sources.setdefault(root_a, set()).add(source)
            continue
        parent[root_b], flip[root_b] = root_a, flip_a ^ flip_b ^ parity
        sources[root_a] = sources.get(root_a, set()) | sources.pop(root_b, set()) | {source}
    groups: dict[str, list[str]] = {}
    for block in preview["blocks"]:
        groups.setdefault(find(block["id"])[0], []).append(block["id"])
    labels = {block["id"]: block["label"] for block in preview["blocks"]}
    membership = {vid: block["id"] for block in preview["blocks"] for vid in block["variant_ids"]}
    linked, info, variants = {}, [], []
    for index, (root, members) in enumerate(((r, m) for r, m in groups.items() if len(m) > 1), 1):
        name = f"linked-{index}"
        for member in members:
            linked[member] = name
        info.append({"phase_set": name, "blocks": [labels[m] for m in members], "sources": sorted(sources.get(root, set()))})
    for variant in locus.variants:
        block = membership.get(variant.id)
        if block not in linked:
            variants.append(variant)
            continue
        orientation = find(block)[1]
        if orientation:
            variant = replace(variant, gt=variant.gt[::-1], alleles=variant.alleles[::-1], cds_alleles=variant.cds_alleles[::-1])
        variants.append(replace(variant, phased=True, phase_set=linked[block]))
    warnings = list(locus.warnings)
    if info:
        merged = sum(len(group["blocks"]) for group in info)
        warnings.append(f"--link-evidence merged {merged} blocks into {len(info)} using clear {' and '.join(sorted({s for g in info for s in g['sources']}))} evidence; the arrangements shown assume it.")
    if skipped:
        warnings.append(f"{skipped} supported relationship{'s' if skipped != 1 else ''} contradicted stronger evidence and {'were' if skipped != 1 else 'was'} not used for linking.")
    return replace(locus, variants=variants, warnings=warnings), info


def label_linked(result: dict, info: list[dict]) -> dict:
    """Name linked blocks after their sources and the blocks they merged."""
    by_set = {group["phase_set"]: group for group in info}
    for block in result.get("blocks", []):
        name = block["label"].removeprefix("Phase set ")
        if name in by_set:
            group = by_set[name]
            block["label"] = f"Linked by {' and '.join(group['sources'])}: " + " + ".join(group["blocks"])
    result["linked_blocks"] = info
    return result
