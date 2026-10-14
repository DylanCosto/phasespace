"""HGVS-style c. and p. labels, including joint edits on one gene copy."""
from pathlib import Path
import shutil
import subprocess

import pytest

from phasespace.analysis import analyze
from phasespace.hgvs import coding_change
from phasespace.inputs import load_locus
from phasespace.models import Evidence, Locus, Transcript, Variant

REFERENCE = "ATGCAACAATGGAAAGGGTAA"  # M Q Q W K G *


def copy_changes(edits, reference=REFERENCE):
    variants = [Variant(name, start + 100, ref, (alt, alt), (1, 1), False, None, start, start + len(ref), ref, (alt, alt), 60, 20)
                for name, start, ref, alt in edits]
    locus = Locus(Transcript("t", "G", "c", "+", ((100, 100 + len(reference)),), reference), "s", variants,
                  [Evidence(100, 100 + len(reference), "reference")])
    return analyze(locus)["scenarios"][0]["haplotypes"][0]["protein_changes"]


@pytest.mark.parametrize("edits,expected", [
    ([("a", 3, "C", "T")], ["p.Gln2Ter"]),
    ([("a", 3, "C", "G")], ["p.Gln2Glu"]),
    ([("a", 3, "C", "G"), ("b", 5, "A", "T")], ["p.Gln2Asp"]),  # Same codon: CAA -> GAT.
    ([("a", 5, "A", "G")], ["p.="]),
    ([("a", 5, "ACAA", "A")], ["p.Gln3del"]),  # Gln Gln -> Gln: the 3' copy is deleted.
    ([("a", 8, "A", "AGGG")], ["p.Gln3_Trp4insGly"]),
    ([("a", 8, "A", "ACAA")], ["p.Gln3dup"]),
    ([("a", 4, "A", "AT")], ["p.Gln2HisfsTer?"]),
    ([("a", 18, "T", "C")], ["p.Ter7GlnextTer?"]),
    ([("a", 0, "A", "G")], ["p.Met1?"]),
    ([("a", 3, "C", "G"), ("b", 12, "A", "T")], ["p.Gln2Glu", "p.Lys5Ter"]),
])
def test_protein_changes_of_one_copy(edits, expected):
    assert copy_changes(edits) == expected


def test_frameshift_counts_codons_to_the_new_stop():
    # ATG AAA CTA AGC GGG TAA (M K L S G *); deleting c.4A reads ATG AAC TAA: Lys2Asn, then a stop.
    assert copy_changes([("a", 2, "GA", "G")], "ATGAAACTAAGCGGGTAA") == ["p.Lys2AsnfsTer2"]
    # ATG GTA AGC TAA (M V S *); deleting c.4G reads ATG TAA: the first changed codon is the stop.
    assert copy_changes([("a", 2, "GG", "G")], "ATGGTAAGCTAA") == ["p.Val2Ter"]


def test_stop_loss_without_a_later_stop_in_the_cds():
    # The CDS ends at its stop codon, so a lost stop reads off the end: extTer?.
    assert copy_changes([("a", 9, "T", "C")], "ATGCAAGGGTAA") == ["p.Ter4GlnextTer?"]


@pytest.mark.parametrize("start,ref,alt,expected", [
    (21, "C", "T", "c.22C>T"), (30, "CCA", "C", "c.32_33del"), (30, "CA", "C", "c.32del"),
    (30, "GTC", "C", "c.31_32del"), (8, "A", "AT", "c.9_10insT"), (8, "A", "TA", "c.8_9insT"), (4, "AC", "GT", "c.5_6delinsGT"),
])
def test_coding_change_labels(start, ref, alt, expected):
    assert coding_change(start, ref, alt) == expected


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_python_and_browser_coding_labels_agree():
    cases = [(21, "C", "T"), (30, "CCA", "C"), (30, "GTC", "C"), (8, "A", "AT"), (8, "A", "TA"), (4, "AC", "GT"), (10, "ACGT", "A"), (3, "G", "GGGA")]
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    script = "const core = require(process.argv[1]); for (const [s, r, a] of JSON.parse(process.argv[2])) console.log(core.changeLabel({cds_start: s, cds_end: s + r.length, cds_ref: r}, a));"
    import json
    out = subprocess.run(["node", "-e", script, str(javascript), json.dumps(cases)], check=True, capture_output=True, text=True).stdout.split()
    assert out == [coding_change(*case) for case in cases]


DATA = Path(__file__).parent / "data" / "hg002_chr20"


def real_changes(vcf, transcript):
    result = analyze(load_locus(DATA / f"{vcf}.HG002.vcf.gz", DATA / "GRCh38.chr20.masked.fa.gz", DATA / "gencode.v46.gtf", transcript), max_blocks=12)
    return [[h["protein_changes"] for h in s["haplotypes"]] for s in result["scenarios"]]


def test_real_protein_changes_on_both_strands():
    assert real_changes("pacbio", "ENST00000371998.8") == [[[], ["p.Gln1255del"]]]  # NCOA3, forward strand.
    assert real_changes("pacbio", "ENST00000372726.5") == [[["p.Asn97del"], ["p.Asn97del"]]]  # TP53TG5, reverse strand.
    assert real_changes("pacbio", "ENST00000382398.4") == [[["p.Pro106ArgfsTer?"]] * 2]  # DEFB126 reads through the stop.
    assert real_changes("trio", "ENST00000217428.7") == [[["p.Leu77Ser"], ["p.Phe54Ser"]]]  # SPINT3, reverse strand.


def test_real_joint_codon_change_in_sirpa():
    # c.298A>G alone is p.Asn100Asp and c.300C>A alone is p.Asn100Lys; trio phase puts both on one copy.
    changes = real_changes("trio", "ENST00000358771.5")
    assert all("p.Asn100Glu" in scenario[1] for scenario in changes)
    assert not any("p.Asn100Asp" in copy or "p.Asn100Lys" in copy for scenario in changes for copy in scenario)
