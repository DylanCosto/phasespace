"""Gene selection, the annotation catalogue and whole-annotation scans."""
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from phasespace.analysis import analyze, summarize
from phasespace.annotation import catalogue, select_transcript
from phasespace.cli import main
from phasespace.inputs import load_locus
from phasespace.models import PhaseSpaceError

DATA = Path(__file__).parent / "data" / "hg002_chr20"
FASTA = DATA / "GRCh38.chr20.masked.fa.gz"
GTF, GFF3 = DATA / "gencode.v46.gtf", DATA / "ensembl112.gff3"
FIXTURE = {"NCOA3": "ENST00000371998.8", "TP53TG5": "ENST00000372726.5", "CDC25B": "ENST00000245960.10",
           "SPINT3": "ENST00000217428.7", "DEFB126": "ENST00000382398.4", "SIRPA": "ENST00000358771.5"}


@pytest.mark.parametrize("annotation", [GTF, GFF3])
@pytest.mark.parametrize("gene", FIXTURE)
def test_gene_selects_its_mane_select_transcript(annotation, gene):
    transcript, reason = select_transcript(annotation, gene)
    assert transcript == FIXTURE[gene]
    assert reason.startswith("MANE Select transcript")


def test_gene_selection_by_gene_id_and_errors():
    assert select_transcript(GTF, "ENSG00000125788")[0] == FIXTURE["DEFB126"]
    with pytest.raises(PhaseSpaceError, match="No protein-coding transcript found"):
        select_transcript(GTF, "NOT_A_GENE")


def test_gene_without_representative_transcript_lists_choices(tmp_path):
    annotation = tmp_path / "plain.gtf"
    annotation.write_text(
        'chr1\tt\ttranscript\t1\t90\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "G"; transcript_type "protein_coding";\n'
        'chr1\tt\ttranscript\t1\t90\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "G"; transcript_type "protein_coding"; tag "Ensembl_canonical";\n')
    assert select_transcript(annotation, "G") == ("T2", "Ensembl canonical transcript T2 selected for gene G.")
    annotation.write_text(annotation.read_text().replace(' tag "Ensembl_canonical";', ""))
    with pytest.raises(PhaseSpaceError, match=r"choose one with --transcript \(T1, T2\)"):
        select_transcript(annotation, "G")


def test_refseq_gff3_mane_select_tag(tmp_path):
    annotation = tmp_path / "refseq.gff3"
    annotation.write_text(
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=gene-G;Name=G\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=rna-NM_2.1;Parent=gene-G;gene=G;transcript_id=NM_2.1\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=rna-NM_1.3;Parent=gene-G;gene=G;transcript_id=NM_1.3;tag=MANE Select\n")
    assert select_transcript(annotation, "G")[0] == "NM_1.3"


@pytest.mark.parametrize("annotation", [GTF, GFF3])
def test_catalogue_assembles_the_same_models_as_single_lookups(annotation):
    models, skipped = catalogue(annotation)
    assert {model.gene: model.id for model in models} == FIXTURE
    assert not skipped
    for model in models:
        single = load_locus(DATA / "pacbio.HG002.vcf.gz", FASTA, GTF, FIXTURE[model.gene]).transcript
        assert model.intervals == single.cds_intervals
        assert model.selection == "MANE Select"


def test_catalogue_reports_unsupported_and_missing_genes(tmp_path):
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(GTF.read_text() + 'chr20\tt\tSelenocysteine\t145700\t145702\t.\t+\t.\tgene_id "ENSG00000125788.6"; transcript_id "ENST00000382398.4";\n')
    models, skipped = catalogue(annotation, {"DEFB126", "SPINT3", "MISSING"})
    assert [model.gene for model in models] == ["SPINT3"]
    assert ("DEFB126", "chr20") == (skipped[0][0], skipped[0][2]) and "Selenocysteine" in skipped[0][1]
    assert ("MISSING", "gene not found", "") in skipped


def scan(tmp_path, *extra, vcf="deepvariant"):
    out = tmp_path / "scan"
    code = main(["scan", "--vcf", str(DATA / f"{vcf}.HG002.vcf.gz"), "--fasta", str(FASTA), "--annotation", str(GTF),
                 "--out-dir", str(out), *extra])
    assert code == 0
    rows = {row["gene"]: row for row in csv.DictReader(open(out / "scan.tsv"), delimiter="\t")}
    return out, rows


@pytest.mark.parametrize("threads", ["1", "2"])
def test_scan_ranks_genes_and_writes_candidate_reports(tmp_path, threads):
    out, rows = scan(tmp_path, "--threads", threads)
    assert set(rows) == set(FIXTURE)
    assert rows["SPINT3"]["tier"] == "1" and rows["SPINT3"]["verdict"] == "mixed"
    assert rows["DEFB126"]["tier"] == "2" and rows["DEFB126"]["disruptive"] == "frameshift,stop_lost"
    assert rows["CDC25B"]["tier"] == "4" and rows["CDC25B"]["changes"] == ""
    assert rows["SIRPA"]["tier"] == "5" and "20 independent phase blocks" in rows["SIRPA"]["status"]
    assert list(rows)[0] == "SPINT3"  # Phase-decides genes come first.
    reports = sorted(path.name for path in (out / "reports").iterdir())
    assert reports == ["DEFB126_ENST00000382398.4.html", "SPINT3_ENST00000217428.7.html"]
    index = (out / "index.html").read_text()
    assert "connect-src 'none'" in index and not re.search(r'(?:src|href)=["\']https?://', index)
    summary = json.loads((out / "scan.json").read_text())
    assert summary["counts"]["phase_decides"] == 1 and summary["counts"]["reports"] == 2


def test_scan_with_trio_phase_and_gene_list(tmp_path):
    _, rows = scan(tmp_path, "--threads", "1", "--max-blocks", "12", "--genes", "SPINT3,SIRPA",
                   "--phase-from", str(DATA / "whatshap.HG002.vcf.gz"))
    assert set(rows) == {"SPINT3", "SIRPA"}
    assert rows["SPINT3"]["tier"] == "3" and rows["SPINT3"]["arrangements"] == "1"  # Resolved: missense on both copies.
    assert rows["SIRPA"]["arrangements"] == "2" and rows["SIRPA"]["phased"] == "19"


def test_scan_gene_file_and_existing_directory(tmp_path):
    genes = tmp_path / "panel.txt"
    genes.write_text("# recessive panel\nNCOA3\n\nDEFB126  # common frameshift\n")
    out, rows = scan(tmp_path, "--threads", "1", "--genes", f"@{genes}", "--reports", "all")
    assert set(rows) == {"NCOA3", "DEFB126"}
    assert len(list((out / "reports").iterdir())) == 2
    assert main(["scan", "--vcf", str(DATA / "pacbio.HG002.vcf.gz"), "--fasta", str(FASTA), "--annotation", str(GTF), "--out-dir", str(out)]) == 2


def test_report_gene_option(tmp_path):
    data = tmp_path / "spint3.json"
    assert main(["report", "--vcf", str(DATA / "pacbio.HG002.vcf.gz"), "--fasta", str(FASTA), "--annotation", str(GFF3),
                 "--gene", "SPINT3", "--out", str(tmp_path / "spint3.html"), "--json", str(data)]) == 0
    result = json.loads(data.read_text())
    assert result["transcript"]["id"] == FIXTURE["SPINT3"]
    assert result["warnings"][0] == "MANE Select transcript ENST00000217428.7 selected for gene SPINT3."


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_python_summary_matches_report_verdict(tmp_path):
    results = []
    for label in ("pacbio", "deepvariant", "trio"):
        for gene, transcript in FIXTURE.items():
            if gene == "SIRPA" and label != "trio":
                continue  # Over the block cap unless phased.
            result = analyze(load_locus(DATA / f"{label}.HG002.vcf.gz", FASTA, GTF, transcript), max_blocks=12)
            results.append({"scenarios": result["scenarios"], "reference": result["transcript"]["reference_cds"], "verdict": summarize(result)["verdict"]})
    data = tmp_path / "results.json"
    data.write_text(json.dumps(results))
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    test = r"""
const assert = require('node:assert/strict');
const core = require(process.argv[1]);
for (const item of JSON.parse(require('node:fs').readFileSync(process.argv[2], 'utf8'))) {
  assert.equal(core.phaseVerdict(item.scenarios, item.reference).tone, item.verdict);
}
"""
    subprocess.run(["node", "-e", test, str(javascript), str(data)], check=True, capture_output=True, text=True)
