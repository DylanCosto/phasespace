"""Read-backed phase evidence: allele typing, fragment counting and real reads."""
import csv
import json
from pathlib import Path
import random

import pysam
import pytest

from phasespace.analysis import analyze
from phasespace.cli import main
from phasespace.inputs import load_locus
from phasespace.reads import add_read_evidence, read_evidence, supported_constraints

DATA = Path(__file__).parent / "data" / "hg002_chr20"
FASTA = DATA / "GRCh38.chr20.masked.fa.gz"
GTF = DATA / "gencode.v46.gtf"

random.seed(5)
SEQUENCE = "".join(random.choice("ACGT") for _ in range(300))
# One CDS [100, 160): codon-aligned, starts ATG, ends TAA, no internal stop.
CDS = "ATG" + "GCT" * 18 + "TAA"
SEQUENCE = SEQUENCE[:100] + CDS + SEQUENCE[160:]
SNV_A, SNV_B, DEL = 109, 130, 145  # 0-based; DEL removes two bases after its anchor.


@pytest.fixture
def locus_files(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + SEQUENCE + "\n")
    pysam.faidx(str(fasta))
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\tt\tCDS\t101\t160\t.\t+\t0\tgene_id "G"; transcript_id "TX"; gene_name "G";\n')
    alt = lambda pos: "A" if SEQUENCE[pos] != "A" else "C"
    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=300>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"chr1\t{SNV_A + 1}\t.\t{SEQUENCE[SNV_A]}\t{alt(SNV_A)}\t.\tPASS\t.\tGT:GQ:DP\t0/1:50:30\n"
        f"chr1\t{SNV_B + 1}\t.\t{SEQUENCE[SNV_B]}\t{alt(SNV_B)}\t.\tPASS\t.\tGT:GQ:DP\t0/1:50:30\n"
        f"chr1\t{DEL + 1}\t.\t{SEQUENCE[DEL:DEL + 4]}\t{SEQUENCE[DEL:DEL + 1] + SEQUENCE[DEL + 3]}\t.\tPASS\t.\tGT:GQ:DP\t0/1:50:30\n"
    )
    return fasta, gtf, vcf, alt


def write_bam(path, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 300}]}
    with pysam.AlignmentFile(str(path) + ".tmp", "wb", header=header) as handle:
        for name, start, sequence, cigar, mapq in reads:
            read = pysam.AlignedSegment(handle.header)
            read.query_name, read.reference_id, read.reference_start = name, 0, start
            read.query_sequence, read.cigarstring, read.mapping_quality = sequence, cigar, mapq
            read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            handle.write(read)
    pysam.sort("-o", str(path), str(path) + ".tmp")
    pysam.index(str(path))
    return path


def haplotype(alt, snv_a: bool, snv_b: bool, deletion: bool, start=90, end=170):
    """A read sequence from one copy, with its CIGAR."""
    bases = list(SEQUENCE[start:end])
    if snv_a:
        bases[SNV_A - start] = alt(SNV_A)
    if snv_b:
        bases[SNV_B - start] = alt(SNV_B)
    if deletion:
        cut = DEL + 1 - start
        del bases[cut:cut + 2]
        return "".join(bases), f"{cut}M2D{len(bases) - cut}M"
    return "".join(bases), f"{len(bases)}M"


def analysed(locus_files):
    fasta, gtf, vcf, _ = locus_files
    return analyze(load_locus(vcf, fasta, gtf, "TX"))


def test_reads_type_snvs_and_deletions_and_count_fragments(locus_files, tmp_path):
    fasta, _, _, alt = locus_files
    result = analysed(locus_files)
    reads = []
    for index in range(9):  # Copy 1: both SNV alternates; copy 2: the deletion.
        reads.append((f"cis{index}", 90, *haplotype(alt, True, True, False), 60))
        reads.append((f"del{index}", 90, *haplotype(alt, False, False, True), 60))
    reads.append(("odd", 90, *haplotype(alt, True, False, False), 60))  # One fragment disagrees.
    reads.append(("lowmapq", 90, *haplotype(alt, False, True, True), 5))  # Ignored.
    bam = write_bam(tmp_path / "reads.bam", reads)
    evidence = read_evidence(result, [bam], fasta)
    ids = [v["id"] for v in result["variants"]]
    assert evidence["fragments"] == 19
    assert [evidence["variants"][vid]["allele_reads"] for vid in ids] == [[9, 10], [10, 9], [10, 9]]
    pairs = {(p["block_a"], p["block_b"]): (p["same"], p["opposite"]) for p in evidence["pairs"]}
    # Allele indexes follow genotype order (0/1). Copy 1 reads are (1, 1, 0) and copy 2 reads (0, 0, 1),
    # so the SNVs share an orientation and the deletion is opposite both; "odd" (1, 0, 0) dissents twice.
    assert pairs == {("block-1", "block-2"): (18, 1), ("block-1", "block-3"): (0, 19), ("block-2", "block-3"): (1, 18)}
    assert {(c["a"], c["b"], c["parity"]) for c in supported_constraints(evidence)} == {("block-1", "block-2", 0), ("block-1", "block-3", 1), ("block-2", "block-3", 1)}


def test_mates_join_into_one_fragment_and_disagreeing_mates_drop_a_site(locus_files, tmp_path):
    fasta, _, _, alt = locus_files
    result = analysed(locus_files)
    first, first_cigar = haplotype(alt, True, False, False, 100, 120)   # Covers SNV_A only.
    second, second_cigar = haplotype(alt, False, True, False, 125, 140)  # Covers SNV_B only.
    clash, clash_cigar = haplotype(alt, False, False, False, 100, 120)   # Same name, other allele at SNV_A.
    bam = write_bam(tmp_path / "pairs.bam", [("pair", 100, first, first_cigar, 60), ("pair", 125, second, second_cigar, 60),
                                             ("bad", 100, first, first_cigar, 60), ("bad", 125, second, second_cigar, 60), ("bad", 100, clash, clash_cigar, 60)])
    evidence = read_evidence(result, [bam], fasta)
    pairs = {(p["block_a"], p["block_b"]): (p["same"], p["opposite"]) for p in evidence["pairs"]}
    assert pairs == {("block-1", "block-2"): (1, 0)}  # Only "pair": its mates show alt at A and alt at B.
    assert evidence["fragments"] == 2  # "bad" still counts, but its clashing site A is dropped.


def test_verdict_with_reads_and_conflicts(locus_files, tmp_path):
    fasta, _, _, alt = locus_files
    result = analysed(locus_files)
    reads = [(f"r{i}", 90, *haplotype(alt, True, True, False), 60) for i in range(5)]
    add_read_evidence(result, [write_bam(tmp_path / "cis.bam", reads)], fasta)
    evidence = result["read_evidence"]
    assert {(c["a"], c["b"], c["parity"]) for c in evidence["supported"]} == {("block-1", "block-2", 0), ("block-1", "block-3", 1), ("block-2", "block-3", 1)}
    assert evidence["compatible_arrangements"] == 1 and evidence["thresholds"] == {"min_count": 3, "min_fraction": 0.9}
    # Below the support threshold nothing is applied.
    few = analysed(locus_files)
    add_read_evidence(few, [write_bam(tmp_path / "two.bam", reads[:2])], fasta)
    assert few["read_evidence"]["supported"] == [] and few["read_evidence"]["compatible_arrangements"] == 4


# Real HG002 reads (GIAB PacBio HiFi Revio and Illumina), checked against trio phasing.

def real(vcf, transcript):
    return analyze(load_locus(DATA / f"{vcf}.HG002.vcf.gz", FASTA, GTF, transcript), max_blocks=12)


@pytest.mark.parametrize("tech,expected", [("hifi", [(0, 53), (52, 0), (0, 51)]), ("illumina", [(63, 0), (0, 41), (0, 34)])])
def test_real_reads_resolve_spint3_as_trio_phasing_does(tech, expected):
    result = real("deepvariant", "ENST00000217428.7")
    add_read_evidence(result, [DATA / f"{tech}.HG002.bam"], FASTA)
    evidence = result["read_evidence"]
    assert [(p["same"], p["opposite"]) for p in evidence["pairs"]] == expected
    assert evidence["compatible_arrangements"] == 1 and evidence["verdict_with"] == "none"
    supported = evidence["supported"]
    chosen = [s for s in result["scenarios"] if all(s["orientations"][c["a"]] ^ s["orientations"][c["b"]] == c["parity"] for c in supported)]
    trio = real("trio", "ENST00000217428.7")
    assert chosen[0]["outcome_key"] == trio["scenarios"][0]["outcome_key"]


def test_long_reads_span_cdc25b_but_short_reads_cannot():
    result = real("deepvariant", "ENST00000245960.10")
    hifi = add_read_evidence(dict(result), [DATA / "hifi.HG002.bam"], FASTA)["read_evidence"]
    assert [(p["same"], p["opposite"]) for p in hifi["pairs"]] == [(28, 0), (19, 0), (15, 0)]
    illumina = add_read_evidence(dict(result), [DATA / "illumina.HG002.bam"], FASTA)["read_evidence"]
    assert illumina["pairs"] == [] and illumina["fragments"] > 100


def test_real_reads_confirm_trio_phase_within_sirpa_and_place_the_unphased_site():
    result = real("trio", "ENST00000358771.5")
    evidence = add_read_evidence(result, [DATA / "hifi.HG002.bam"], FASTA)["read_evidence"]
    assert evidence["within_blocks"] == {"agree": 6969, "disagree": 0}
    assert [(p["same"], p["opposite"]) for p in evidence["pairs"]] == [(40, 2)]
    assert evidence["compatible_arrangements"] == 1


def test_cli_reads_in_report_and_scan(tmp_path):
    data = tmp_path / "spint3.json"
    common = ["--vcf", str(DATA / "deepvariant.HG002.vcf.gz"), "--fasta", str(FASTA), "--annotation", str(GTF), "--reads", str(DATA / "hifi.HG002.bam")]
    assert main(["report", *common, "--gene", "SPINT3", "--out", str(tmp_path / "spint3.html"), "--json", str(data)]) == 0
    assert json.loads(data.read_text())["read_evidence"]["verdict_with"] == "none"
    assert main(["scan", *common, "--genes", "SPINT3,CDC25B", "--threads", "1", "--out-dir", str(tmp_path / "scan")]) == 0
    rows = {row["gene"]: row for row in csv.DictReader(open(tmp_path / "scan" / "scan.tsv"), delimiter="\t")}
    assert (rows["SPINT3"]["verdict"], rows["SPINT3"]["verdict_with_reads"], rows["SPINT3"]["read_supported_relationships"]) == ("mixed", "none", "3")
    assert rows["CDC25B"]["verdict_with_reads"] == "good"
