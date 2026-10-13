"""Regression checks on real GIAB HG002 calls (see tests/data/hg002_chr20/README.md).

Each case was cross-checked against bcftools csq --local-csq on the same
transcripts; per-variant consequences agreed except that a frameshift reading
through the stop codon is also labelled stop_lost here.
"""
import json
from pathlib import Path
import re

import pytest

from phasespace.analysis import analyze
from phasespace.cli import main
from phasespace.inputs import load_locus

DATA = Path(__file__).parent / "data" / "hg002_chr20"
FASTA = DATA / "GRCh38.chr20.masked.fa.gz"
INPUTS = ["pacbio", "deepvariant", "giab"]
ANNOTATIONS = [("gencode.v46.gtf", lambda tx: tx), ("ensembl112.gff3", lambda tx: tx.split(".")[0])]

# transcript: gene, strand, CDS length incl. stop, genotypes, blocks, arrangements,
# protein outcomes, distinct per-copy consequence labels, reference first stop
EXPECTED = {
    "ENST00000371998.8": ("NCOA3", "+", 4275, [(0, 1)], 1, 1, 1, {("inframe_deletion",), ("reference",)}, 1425),
    "ENST00000372726.5": ("TP53TG5", "-", 873, [(1, 1)], 0, 1, 1, {("inframe_deletion",)}, 291),
    "ENST00000245960.10": ("CDC25B", "+", 1743, [(0, 1)] * 3, 3, 4, 1, {("reference",), ("synonymous",)}, 581),
    "ENST00000217428.7": ("SPINT3", "-", 270, [(0, 1)] * 3, 3, 4, 2, {("missense",), ("reference",), ("synonymous",)}, 90),
    "ENST00000382398.4": ("DEFB126", "+", 336, [(1, 1)], 0, 1, 1, {("frameshift", "stop_lost")}, 112),
}
# Coding bases without usable evidence; the variant-only GIAB VCF leaves almost all unknown.
UNKNOWN = {
    "pacbio": {"NCOA3": 2, "TP53TG5": 3, "CDC25B": 4, "SPINT3": 0, "DEFB126": 0},
    "deepvariant": {"NCOA3": 5, "TP53TG5": 2, "CDC25B": 0, "SPINT3": 1, "DEFB126": 0},
}


def load(label, annotation, transcript):
    name, to_id = annotation
    return load_locus(DATA / f"{label}.HG002.vcf.gz", FASTA, DATA / name, to_id(transcript))


@pytest.mark.parametrize("annotation", ANNOTATIONS, ids=[name for name, _ in ANNOTATIONS])
@pytest.mark.parametrize("label", INPUTS)
@pytest.mark.parametrize("transcript", EXPECTED)
def test_real_calls_give_stable_results_across_callers_and_annotations(transcript, label, annotation):
    gene, strand, length, genotypes, blocks, arrangements, proteins, labels, ref_stop = EXPECTED[transcript]
    locus = load(label, annotation, transcript)
    result = analyze(locus)
    assert (locus.transcript.gene, locus.transcript.strand, len(locus.transcript.reference_cds)) == (gene, strand, length)
    # GTF stop_codon features are joined to the CDS, so both formats end in a stop.
    assert locus.transcript.reference_cds[-3:] in ("TAA", "TAG", "TGA")
    assert not any("terminal stop" in warning for warning in result["warnings"])
    assert [v.gt for v in locus.variants] == genotypes
    assert (len(result["blocks"]), result["scenario_count"], result["protein_outcome_count"]) == (blocks, arrangements, proteins)
    assert {tuple(h["consequences"]) for s in result["scenarios"] for h in s["haplotypes"]} == labels
    assert result["reference_first_stop"] == ref_stop
    unknown = sum(e.end - e.start for e in locus.evidence if e.status == "unknown")
    if label == "giab":
        assert unknown >= length - 4  # Variant-only VCF: absent records never become reference.
    else:
        assert unknown == UNKNOWN[label][gene]


@pytest.mark.parametrize("transcript", EXPECTED)
def test_gtf_and_gff3_describe_the_same_coding_sequence(transcript):
    gtf = load("pacbio", ANNOTATIONS[0], transcript).transcript
    gff3 = load("pacbio", ANNOTATIONS[1], transcript).transcript
    assert gtf.cds_intervals == gff3.cds_intervals
    assert gtf.reference_cds == gff3.reference_cds


def test_versioned_and_unversioned_ids_resolve_in_both_formats():
    for name, identifier in [("gencode.v46.gtf", "ENST00000371998"), ("ensembl112.gff3", "ENST00000371998.8"), ("ensembl112.gff3", "transcript:ENST00000371998")]:
        locus = load_locus(DATA / "pacbio.HG002.vcf.gz", FASTA, DATA / name, identifier)
        assert locus.transcript.gene == "NCOA3"


def test_deepvariant_reference_block_abutting_a_heterozygous_deletion(tmp_path):
    # DeepVariant starts the next reference block after the deleted bases.
    locus = load("deepvariant", ANNOTATIONS[0], "ENST00000371998.8")
    deletion = locus.variants[0]
    assert (deletion.ref, deletion.alleles) == ("ACAG", ("ACAG", "A"))
    deleted = {base for e in locus.evidence if e.status == "variant" for base in range(e.start, e.end)}
    assert set(range(deletion.pos, deletion.pos + 4)) <= deleted


def test_cli_report_on_real_data(tmp_path):
    out, data = tmp_path / "spint3.html", tmp_path / "spint3.json"
    code = main(["report", "--vcf", str(DATA / "pacbio.HG002.vcf.gz"), "--fasta", str(FASTA),
                 "--annotation", str(DATA / "gencode.v46.gtf"), "--transcript", "ENST00000217428.7",
                 "--out", str(out), "--json", str(data)])
    assert code == 0
    result = json.loads(data.read_text())
    assert result["sample"] == "HG002" and result["scenario_count"] == 4
    embedded = re.search(r'<script id="phasespace-data" type="application/json">(.*?)</script>', out.read_text(), re.S)
    assert json.loads(embedded.group(1))["report_encoding"] == "edits"
