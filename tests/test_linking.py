"""--link-evidence: merging blocks that reads or parents clearly link."""
import csv
from pathlib import Path

import pytest

from phasespace.analysis import analyze
from phasespace.cli import main
from phasespace.inputs import load_locus
from phasespace.linking import label_linked, link_blocks
from phasespace.models import PhaseSpaceError
from test_reads import SNV_A, SNV_B, haplotype, locus_files, write_bam  # noqa: F401 (fixture)

DATA = Path(__file__).parent / "data" / "hg002_chr20"
FASTA = DATA / "GRCh38.chr20.masked.fa.gz"
GTF = DATA / "gencode.v46.gtf"
SIRPA = "ENST00000358771.5"


def sirpa():
    return load_locus(DATA / "deepvariant.HG002.vcf.gz", FASTA, GTF, SIRPA)


def test_linking_brings_sirpa_under_the_cap_and_sources_agree():
    with pytest.raises(PhaseSpaceError, match="20 independent phase blocks"):
        analyze(sirpa(), max_blocks=12)
    outcomes = {}
    for name, options in {"hifi": {"reads": [DATA / "hifi.HG002.bam"], "fasta": FASTA},
                          "illumina": {"reads": [DATA / "illumina.HG002.bam"], "fasta": FASTA},
                          "parents": {"parents": ("HG003", "HG004"), "parents_from": DATA / "whatshap.HG002.vcf.gz"}}.items():
        linked, info = link_blocks(sirpa(), **options)
        assert [len(group["blocks"]) for group in info] == [20]
        result = label_linked(analyze(linked), info)
        assert result["scenario_count"] == 1 and result["summary"]["verdict"] == "good"
        assert result["blocks"][0]["label"].startswith("Linked by ")
        assert any("--link-evidence merged 20 blocks into 1" in w for w in result["warnings"])
        outcomes[name] = result["scenarios"][0]["outcome_key"]
    assert len(set(outcomes.values())) == 1  # HiFi, Illumina and parents agree on the same pair of copies.


def test_linking_skips_a_weaker_contradicting_relationship(locus_files, tmp_path):
    fasta, gtf, vcf, alt = locus_files
    locus = load_locus(vcf, fasta, gtf, "TX")
    reads = []
    for i in range(10):  # A and B together (alternate alleles on one copy).
        reads.append((f"ab{i}", 100, *haplotype(alt, True, True, False, 100, 135), 60))
    for i in range(10):  # B and the deletion together.
        reads.append((f"bd{i}", 125, *haplotype(alt, False, True, True, 125, 165), 60))
    for i in range(4):  # Mate pairs: A alternate with the reference at the deletion -- contradicts the two above.
        first, cigar = haplotype(alt, True, False, False, 100, 120)
        second, cigar2 = haplotype(alt, False, False, False, 140, 165)
        reads += [(f"mate{i}", 100, first, cigar, 60), (f"mate{i}", 140, second, cigar2, 60)]
    linked, info = link_blocks(locus, reads=[write_bam(tmp_path / "mixed.bam", reads)], fasta=fasta)
    assert [len(group["blocks"]) for group in info] == [3]
    assert any("1 supported relationship contradicted stronger evidence" in w for w in linked.warnings)
    # Strongest-first: A-B same and B-deletion same, so the deletion joins A's copy.
    result = analyze(linked)
    copies = [h["variant_ids"] for h in result["scenarios"][0]["haplotypes"]]
    assert sorted(len(c) for c in copies) == [0, 3]


def test_cli_link_evidence_in_report_and_scan(tmp_path):
    common = ["--vcf", str(DATA / "deepvariant.HG002.vcf.gz"), "--fasta", str(FASTA), "--annotation", str(GTF), "--max-blocks", "12"]
    assert main(["report", *common, "--gene", "SIRPA", "--link-evidence", "--out", str(tmp_path / "x.html")]) == 2  # Needs evidence.
    assert main(["report", *common, "--gene", "SIRPA", "--link-evidence", "--reads", str(DATA / "hifi.HG002.bam"), "--out", str(tmp_path / "sirpa.html")]) == 0
    assert main(["scan", *common, "--genes", "SIRPA,SPINT3", "--threads", "1", "--link-evidence", "--reads", str(DATA / "hifi.HG002.bam"),
                 "--out-dir", str(tmp_path / "scan")]) == 0
    rows = {row["gene"]: row for row in csv.DictReader(open(tmp_path / "scan" / "scan.tsv"), delimiter="\t")}
    assert (rows["SIRPA"]["tier"], rows["SIRPA"]["linked_blocks"], rows["SIRPA"]["arrangements"]) == ("4", "20", "1")
    assert (rows["SPINT3"]["linked_blocks"], rows["SPINT3"]["verdict"]) == ("3", "none")
