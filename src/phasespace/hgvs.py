"""HGVS-style labels for coding (c.) and protein (p.) changes.

Substitutions follow HGVS exactly. Indels keep the VCF representation and are
not shifted 3'. Protein changes describe each gene copy as a whole, from its
jointly edited coding sequence: two edits in one codon give the amino acid
actually produced, and nothing after the first stop is described.
"""
from __future__ import annotations

THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa",
}


def _aa(amino: str) -> str:
    return THREE_LETTER.get(amino, "Xaa")


def _aas(aminos: str) -> str:
    return "".join(_aa(amino) for amino in aminos)


def coding_change(cds_start: int, cds_ref: str, allele: str) -> str:
    """c. label for one allele (mirrors report.js changeLabel)."""
    ref, alt, start = cds_ref.upper(), allele.upper(), cds_start
    span = lambda a, b: str(a) if a == b else f"{a}_{b}"
    if len(ref) == len(alt):
        return f"c.{start + 1}{ref}>{alt}" if len(ref) == 1 else f"c.{span(start + 1, start + len(ref))}delins{alt}"
    if len(alt) < len(ref) and ref.startswith(alt):
        return f"c.{span(start + len(alt) + 1, start + len(ref))}del"
    if len(alt) < len(ref) and ref.endswith(alt):
        return f"c.{span(start + 1, start + len(ref) - len(alt))}del"
    if len(alt) > len(ref) and alt.startswith(ref):
        return f"c.{start + len(ref)}_{start + len(ref) + 1}ins{alt[len(ref):]}"
    if len(alt) > len(ref) and alt.endswith(ref):
        return f"c.{start}_{start + 1}ins{alt[:len(alt) - len(ref)]}"
    return f"c.{span(start + 1, start + len(ref))}delins{alt}"


def protein_changes(ref_full: str, ref_first_stop: int | None, full: str, first_stop: int | None,
                    alt_to_ref: dict[int, int], frameshift: bool) -> list[str]:
    """p. changes of one gene copy relative to the reference protein.

    ``alt_to_ref`` maps alternate codon indexes to the reference codon they
    descend from intact; unmapped codons are inserted or frameshifted.
    Returns ["p.="] when the protein is unchanged.
    """
    reference = ref_full[:ref_first_stop] if ref_first_stop else ref_full
    alternate = full[:first_stop] if first_stop else full
    if reference[:1] == "M" and alternate[:1] != "M":
        return ["p.Met1?"]
    if frameshift:
        return _frameshift(reference, alternate, first_stop, alt_to_ref)
    changes: list[str] = []
    k = r = 0
    while k < len(alternate) or r < len(reference):
        if k < len(alternate) and alt_to_ref.get(k) == r:
            ref_aa, alt_aa = reference[r], alternate[k]
            if ref_aa == "*" and alt_aa != "*":
                tail = f"Ter{first_stop - 1 - k}" if first_stop else "Ter?"
                changes.append(f"p.Ter{r + 1}{_aa(alt_aa)}ext{tail}")
                break
            if ref_aa != alt_aa:
                changes.append(f"p.{_aa(ref_aa)}{r + 1}{_aa(alt_aa)}")
            if alt_aa == "*":
                break
            k, r = k + 1, r + 1
            continue
        # An unaligned stretch: find the next codon that descends intact from the reference.
        k2 = next((j for j in range(k, len(alternate)) if alt_to_ref.get(j) is not None and alt_to_ref[j] >= r), len(alternate))
        r2 = alt_to_ref[k2] if k2 < len(alternate) else len(reference)
        inserted, deleted = alternate[k:k2], reference[r:r2]
        change = _indel(reference, r, deleted, inserted)
        if change:
            changes.append(change)
        if "*" in inserted:
            break
        k, r = k2, r2
    return changes or ["p.="]


def _indel(reference: str, start: int, deleted: str, inserted: str) -> str | None:
    """Name an unaligned stretch after trimming residues it shares at both ends.

    Trimming the shared start first places deletions as far 3' as this
    stretch allows; repeats beyond the stretch are not shifted further.
    """
    shared = 0
    while shared < min(len(deleted), len(inserted)) and deleted[shared] == inserted[shared]:
        shared += 1
    deleted, inserted, start = deleted[shared:], inserted[shared:], start + shared
    shared = 0
    while shared < min(len(deleted), len(inserted)) and deleted[-1 - shared] == inserted[-1 - shared]:
        shared += 1
    if shared:
        deleted, inserted = deleted[:-shared], inserted[:-shared]
    end = start + len(deleted)
    position = lambda index: f"{_aa(reference[index])}{index + 1}"
    span = lambda a, b: position(a) if b - a <= 1 else f"{position(a)}_{position(b - 1)}"
    if not deleted and not inserted:
        return None
    if len(deleted) == 1 and len(inserted) == 1:
        return f"p.{position(start)}{_aa(inserted)}"
    if deleted and not inserted:
        return f"p.{span(start, end)}del"
    if inserted and not deleted:
        if start >= len(inserted) and reference[start - len(inserted):start] == inserted:
            return f"p.{span(start - len(inserted), start)}dup"
        if 0 < start < len(reference):
            return f"p.{position(start - 1)}_{position(start)}ins{_aas(inserted)}"
        return f"p.{span(max(start - 1, 0), max(start, 1))}delins{_aas(inserted)}"
    return f"p.{span(start, end)}delins{_aas(inserted)}"


def _frameshift(reference: str, alternate: str, first_stop: int | None, alt_to_ref: dict[int, int]) -> list[str]:
    """Describe a frameshift by its first changed amino acid and the new stop."""
    k = 0
    while k < len(alternate) and alt_to_ref.get(k) == k and k < len(reference) and alternate[k] == reference[k]:
        k += 1
    # Continue past amino acids the shifted frame happens to keep.
    r = alt_to_ref.get(k - 1, k - 1) + 1 if k else 0
    while k < len(alternate) and r < len(reference) and alternate[k] == reference[r]:
        k, r = k + 1, r + 1
    if k >= len(alternate) or r >= len(reference):
        return ["p.="] if alternate == reference else [f"p.{_aa(reference[min(r, len(reference) - 1)])}{min(r, len(reference) - 1) + 1}fs"]
    if alternate[k] == "*":
        return [f"p.{_aa(reference[r])}{r + 1}Ter"]
    tail = f"Ter{first_stop - k}" if first_stop else "Ter?"
    return [f"p.{_aa(reference[r])}{r + 1}{_aa(alternate[k])}fs{tail}"]
