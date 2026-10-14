"""Transmission (parent-of-origin) evidence."""
import csv
import json
from pathlib import Path

import pysam
import pytest

from phasespace.analysis import analyze
from phasespace.cli import main
from phasespace.inheritance import _paternal_index, add_inheritance_evidence
from phasespace.inputs import load_locus
from phasespace.models import PhaseSpaceError

DATA = Path(__file__).parent / "data" / "hg002_chr20"
FASTA = DATA / "GRCh38.chr20.masked.fa.gz"
GTF = DATA / "gencode.v46.gtf"


@pytest.mark.parametrize("father,mother,expected", [
    (("A", "G"), ("A", "A"), (1, "informative")),    # Alt from the father.
    (("A", "A"), ("A", "G"), (0, "informative")),    # Alt from the mother, so ref from the father.
    (("G", "G"), ("A", "A"), (1, "informative")),
    (("A", "G"), ("A", "G"), (None, "uninformative")),  # Both heterozygous.
    (("A", "G"), None, (None, "uninformative")),     # Mother unknown: either way fits.
    (("A", "A"), None, (0, "informative")),          # Father can only give ref.
    (("A", "A"), ("A", "A"), (None, "conflict")),    # Neither parent carries G.
    (None, None, (None, "unknown")),
])
def test_paternal_allele_from_parent_genotypes(father, mother, expected):
    assert _paternal_index(("A", "G"), father, mother) == expected


TRIO_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=90>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Min DP">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHILD\tDAD\tMOM
"""


@pytest.fixture
def trio(tmp_path):
    sequence = list("A" * 90)
    sequence[12:18] = "ATGAAA"
    sequence[30:39] = "CCCGGGTAA"
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(reference))
    annotation = tmp_path / "genes.gtf"
    annotation.write_text('chr1\tt\tCDS\t13\t18\t.\t+\t0\tgene_id "G"; transcript_id "TX";\nchr1\tt\tCDS\t31\t39\t.\t+\t0\tgene_id "G"; transcript_id "TX";\n')

    def write(rows):
        path = tmp_path / "trio.vcf"
        path.write_text(TRIO_HEADER + "".join(row + "\n" for row in rows))
        pysam.tabix_compress(str(path), str(path) + ".gz", force=True)
        pysam.tabix_index(str(path) + ".gz", preset="vcf", force=True)
        return str(path) + ".gz"

    return reference, annotation, write


def site(pos, ref, alt, child, dad, mom, fmt="GT:GQ:DP"):
    return f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t{child}\t{dad}\t{mom}"


def test_trans_compound_heterozygote_from_two_parents(trio):
    reference, annotation, write = trio
    path = write([site(16, "A", "T", "0/1:50:20", "0/1:50:20", "0/0:50:20"), site(32, "C", "T", "0/1:50:20", "0/0:50:20", "0/1:50:20")])
    result = analyze(load_locus(path, reference, annotation, "TX", sample="CHILD"))
    evidence = add_inheritance_evidence(result, path, "DAD", "MOM")["inheritance_evidence"]
    assert [evidence["variants"][v["id"]]["paternal_allele"] for v in result["variants"]] == [1, 0]
    assert [(p["same"], p["opposite"]) for p in evidence["pairs"]] == [(0, 1)]  # Alt alleles on different copies.
    assert evidence["supported"] == [{"a": "block-1", "b": "block-2", "parity": 1}]
    assert evidence["compatible_arrangements"] == 1


def test_parent_reference_block_conflicts_and_missing_depth(trio):
    reference, annotation, write = trio
    rows = [
        "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=20\tGT:GQ:MIN_DP\t./.:.:.\t0/0:50:12\t./.:.:.",  # Dad hom-ref block over c.4.
        site(16, "A", "T", "0/1:50:20", "./.:.:.", "0/1:50:."),  # Mom het, depth missing (as in merged trios).
        site(32, "C", "T", "0/1:50:20", "0/0:50:20", "0/0:50:20"),  # Neither parent carries T.
    ]
    path = write(rows)
    result = analyze(load_locus(path, reference, annotation, "TX", sample="CHILD"))
    evidence = add_inheritance_evidence(result, path, "DAD", "MOM")["inheritance_evidence"]
    statuses = {vid: (v["status"], v["paternal_allele"], v["father"]) for vid, v in evidence["variants"].items()}
    assert statuses == {"chr1:16:A>T": ("informative", 0, "A/A"), "chr1:32:C>T": ("conflict", None, "C/C")}
    assert evidence["conflicts"] == 1 and evidence["pairs"] == []


def test_missing_parent_sample_is_an_error(trio):
    reference, annotation, write = trio
    path = write([site(16, "A", "T", "0/1:50:20", "0/1:50:20", "0/0:50:20")])
    result = analyze(load_locus(path, reference, annotation, "TX", sample="CHILD"))
    with pytest.raises(PhaseSpaceError, match="Parent sample NOBODY not in"):
        add_inheritance_evidence(result, path, "NOBODY", "MOM")


# Real HG002 trio: parents HG003 (father) and HG004 (mother) from the merged trio VCF.

def test_real_inheritance_resolves_spint3_like_trio_phasing():
    result = analyze(load_locus(DATA / "deepvariant.HG002.vcf.gz", FASTA, GTF, "ENST00000217428.7"))
    evidence = add_inheritance_evidence(result, DATA / "whatshap.HG002.vcf.gz", "HG003", "HG004")["inheritance_evidence"]
    assert sorted(v["status"] for v in evidence["variants"].values()) == ["informative"] * 3
    assert evidence["compatible_arrangements"] == 1 and evidence["verdict_with"] == "none"
    trio = analyze(load_locus(DATA / "trio.HG002.vcf.gz", FASTA, GTF, "ENST00000217428.7"))
    chosen = [s for s in result["scenarios"] if all(s["orientations"][c["a"]] ^ s["orientations"][c["b"]] == c["parity"] for c in evidence["supported"])]
    assert chosen[0]["outcome_key"] == trio["scenarios"][0]["outcome_key"]


def test_cli_parents_in_report_and_scan(tmp_path):
    common = ["--vcf", str(DATA / "deepvariant.HG002.vcf.gz"), "--fasta", str(FASTA), "--annotation", str(GTF),
              "--parents", "HG003,HG004", "--parents-from", str(DATA / "whatshap.HG002.vcf.gz")]
    data = tmp_path / "spint3.json"
    assert main(["report", *common, "--gene", "SPINT3", "--out", str(tmp_path / "spint3.html"), "--json", str(data)]) == 0
    assert json.loads(data.read_text())["inheritance_evidence"]["verdict_with"] == "none"
    assert main(["scan", *common, "--genes", "SPINT3,CDC25B", "--threads", "1", "--out-dir", str(tmp_path / "scan")]) == 0
    rows = {row["gene"]: row for row in csv.DictReader(open(tmp_path / "scan" / "scan.tsv"), delimiter="\t")}
    assert (rows["SPINT3"]["verdict_with_parents"], rows["SPINT3"]["mendelian_conflicts"]) == ("none", "0")
    assert rows["CDC25B"]["verdict_with_parents"] == "good"
    assert main(["report", *common[:6], "--parents", "HG003", "--gene", "SPINT3", "--out", str(tmp_path / "x.html")]) == 2
