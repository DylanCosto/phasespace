import json

from phasespace.cli import main
from phasespace.demo import write_demo


def test_demo_runs_end_to_end_with_expected_stops(tmp_path):
    html, result = tmp_path / "demo.html", tmp_path / "demo.json"
    assert main(["demo", "--out", str(html), "--json", str(result)]) == 0
    data = json.loads(result.read_text())
    assert data["synthetic"] is True
    assert data["scenario_count"] == 2
    assert data["outcome_count"] == 2
    assert sorted(tuple(sorted(h["first_stop"] for h in s["haplotypes"])) for s in data["scenarios"]) == [(8, 21), (8, 32)]
    assert sum(e["end"] - e["start"] for e in data["evidence"] if e["status"] == "unknown") == 6
    assert "<html" in html.read_text()


def test_plain_vcf_gaps_are_unknown(tmp_path):
    files = write_demo(tmp_path / "inputs")
    result = tmp_path / "plain.json"
    assert main(["report", "--vcf", str(files["plain_vcf"]), "--fasta", str(files["fasta"]), "--annotation", str(files["annotation"]), "--transcript", "DEMO1.1", "--out", str(tmp_path / "plain.html"), "--json", str(result)]) == 0
    data = json.loads(result.read_text())
    assert sum(e["end"] - e["start"] for e in data["evidence"] if e["status"] == "unknown") == 94
    assert not any(e["status"] == "reference" for e in data["evidence"])


def test_existing_output_requires_force(tmp_path):
    output = tmp_path / "existing.html"
    output.write_text("keep me")
    assert main(["demo", "--out", str(output)]) == 2
    assert output.read_text() == "keep me"


def test_identical_output_paths_rejected(tmp_path):
    output = tmp_path / "both.html"
    assert main(["demo", "--out", str(output), "--json", str(output)]) == 2
    assert not output.exists()


def test_nan_quality_threshold_rejected(tmp_path):
    assert main(["report", "--vcf", "none", "--fasta", "none", "--annotation", "none", "--transcript", "none", "--min-gq", "nan", "--out", str(tmp_path / "out")]) == 2


def test_never_overwrite_input_even_with_force(tmp_path):
    files = write_demo(tmp_path / "inputs")
    original = files["fasta"].read_bytes()
    assert main(["report", "--vcf", str(files["vcf"]), "--fasta", str(files["fasta"]), "--annotation", str(files["annotation"]), "--transcript", "DEMO1.1", "--out", str(files["fasta"]), "--force"]) == 2
    assert files["fasta"].read_bytes() == original
