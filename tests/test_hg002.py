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
from phasespace.models import PhaseSpaceError

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
    "pacbio": {"NCOA3": 2, "TP53TG5": 2, "CDC25B": 4, "SPINT3": 0, "DEFB126": 0},
    "deepvariant": {"NCOA3": 0, "TP53TG5": 2, "CDC25B": 0, "SPINT3": 0, "DEFB126": 0},
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


# Trio phasing: the `trio` input is the `deepvariant` gVCF with phase from transmission.
GTF = ANNOTATIONS[0]
KEEPS_PROTEIN = {"reference", "synonymous", "protein_sequence_unchanged"}


@pytest.mark.parametrize("transcript", EXPECTED)
def test_trio_phase_only_removes_arrangements(transcript):
    unphased = analyze(load("deepvariant", GTF, transcript))
    phased_locus = load("trio", GTF, transcript)
    phased = analyze(phased_locus)
    # Phase can only exclude arrangements; it never creates a new outcome.
    assert {s["outcome_key"] for s in phased["scenarios"]} <= {s["outcome_key"] for s in unphased["scenarios"]}
    assert phased["scenario_count"] <= unphased["scenario_count"]
    assert all(v.phase_set == "157484" for v in phased_locus.variants if v.phased)


def test_trio_phase_resolves_spint3_to_missense_on_both_copies():
    unphased = analyze(load("deepvariant", GTF, "ENST00000217428.7"))
    phased = analyze(load("trio", GTF, "ENST00000217428.7"))
    keeps = [any(set(h["consequences"]) <= KEEPS_PROTEIN for h in s["haplotypes"]) for s in unphased["scenarios"]]
    assert True in keeps and False in keeps  # Unphased, phase decides.
    assert [b["label"] for b in phased["blocks"]] == ["Phase set 157484"]
    assert phased["scenario_count"] == 1
    assert [h["consequences"] for h in phased["scenarios"][0]["haplotypes"]] == [["missense"], ["missense"]]


def test_trio_phase_brings_sirpa_under_the_block_cap():
    with pytest.raises(PhaseSpaceError, match="20 independent phase blocks"):
        analyze(load("deepvariant", GTF, "ENST00000358771.5"))
    result = analyze(load("trio", GTF, "ENST00000358771.5"))
    assert [len(b["variant_ids"]) for b in result["blocks"]] == [19, 1]
    assert result["scenario_count"] == 2
    # One copy carries the 19 phased alternate alleles, including an in-frame deletion.
    for scenario in result["scenarios"]:
        assert "inframe_deletion" in scenario["haplotypes"][1]["consequences"]


@pytest.mark.parametrize("transcript", [*EXPECTED, "ENST00000358771.5"])
def test_phase_from_matches_a_gvcf_with_phase_written_in(transcript):
    # The WhatsHap trio VCF numbers alleles across all three samples; matching is by sequence.
    direct = load_locus(DATA / "deepvariant.HG002.vcf.gz", FASTA, DATA / "gencode.v46.gtf", transcript,
                        phase_from=DATA / "whatshap.HG002.vcf.gz")
    written = load("trio", GTF, transcript)
    summary = lambda locus: [(v.id, v.gt, v.cds_alleles, v.phased, v.phase_set) for v in locus.variants]
    assert summary(direct) == summary(written)
    assert direct.evidence == written.evidence
    assert analyze(direct, max_blocks=12)["scenarios"] == analyze(written, max_blocks=12)["scenarios"]


def test_cli_phase_from_on_real_data(tmp_path):
    out, data = tmp_path / "sirpa.html", tmp_path / "sirpa.json"
    args = ["report", "--vcf", str(DATA / "deepvariant.HG002.vcf.gz"), "--fasta", str(FASTA),
            "--annotation", str(DATA / "gencode.v46.gtf"), "--transcript", "ENST00000358771.5",
            "--out", str(out), "--json", str(data), "--max-blocks", "12"]
    assert main(args) == 2  # Unphased: 20 independent blocks exceed the cap.
    assert main(args + ["--phase-from", str(DATA / "whatshap.HG002.vcf.gz")]) == 0
    result = json.loads(data.read_text())
    assert result["scenario_count"] == 2
    assert result["source"] == "deepvariant.HG002.vcf.gz + phase from whatshap.HG002.vcf.gz"
    assert any("19 of 20 heterozygous coding calls phased" in w for w in result["warnings"])


def test_deepvariant_refcall_counts_as_reference_evidence():
    # Five NCOA3 coding bases are covered only by DeepVariant RefCall records (0/0, GQ >= 20).
    locus = load("deepvariant", GTF, "ENST00000371998.8")
    refcall = [e for e in locus.evidence if "RefCall" in e.detail]
    assert sum(e.end - e.start for e in refcall) == 5
    assert all(e.status == "reference" for e in refcall)


def test_ensembl_contig_names_without_chr_prefix(tmp_path):
    # Ensembl names the contig "20"; the reference and VCF use "chr20".
    plain = tmp_path / "ensembl.gff3"
    plain.write_text("".join(line[3:] if line.startswith("chr20\t") else line for line in (DATA / "ensembl112.gff3").read_text().splitlines(True)))
    locus = load_locus(DATA / "pacbio.HG002.vcf.gz", FASTA, plain, "ENST00000371998")
    assert locus.transcript.contig == "chr20"
    assert locus.variants[0].id == "chr20:47651092:ACAG>A"
    assert any("Contig names differ" in w for w in locus.warnings)
    assert locus.evidence == load("pacbio", GTF, "ENST00000371998.8").evidence
