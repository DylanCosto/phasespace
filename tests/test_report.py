"""Report safety and the exact browser-side filtering semantics."""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from phasespace.report import _script_json, write_report


def minimal_result():
    return {
        "schema_version": 1,
        "sample": "sample 1",
        "transcript": {
            "id": "TX1", "gene": "GENE1", "contig": "chr1", "strand": "+",
            "cds_intervals": [[0, 9]], "reference_cds": "ATGAAATAA",
        },
        "variants": [], "evidence": [], "warnings": [], "blocks": [],
        "scenarios": [{"id": "S1", "orientations": {}, "haplotypes": [
            {"variant_ids": [], "cds": "ATGAAATAA", "protein": "MK", "protein_length": 2,
             "first_stop": 3, "consequences": []},
            {"variant_ids": [], "cds": "ATGAAATAA", "protein": "MK", "protein_length": 2,
             "first_stop": 3, "consequences": []},
        ], "outcome_key": "same", "unaffected_copies": 2}],
        "scenario_count": 1, "outcome_count": 1, "reference_protein": "MK",
        "reference_protein_length": 2, "reference_first_stop": 3,
        "informative_pairs": [], "source": "local.vcf",
    }


def test_report_is_portable_and_preserves_analysis(tmp_path):
    result = minimal_result()
    output = tmp_path / "report.html"
    write_report(result, output)
    html = output.read_text()
    match = re.search(r'<script id="phasespace-data" type="application/json">(.*?)</script>', html, re.S)
    assert match and json.loads(match.group(1)) == result
    assert "PHASESPACE_DATA" not in html
    assert "PHASESPACE_CSS" not in html
    assert "PHASESPACE_JS" not in html
    assert not re.search(r'<(?:script|link)\b[^>]*(?:src|href)=["\']https?://', html)
    assert "connect-src 'none'" in html
    assert "reference-identical" in html.lower()


def test_metadata_cannot_escape_json_or_become_markup(tmp_path):
    result = minimal_result()
    hostile = '</script><img src=x onerror="alert(1)">&\u2028\u2029'
    result["sample"] = hostile
    result["source"] = "/* PHASESPACE_JS */ PHASESPACE_DATA"
    result["warnings"] = [hostile]
    output = tmp_path / "hostile.html"
    write_report(result, output)
    html = output.read_text()
    assert hostile not in html
    match = re.search(r'<script id="phasespace-data" type="application/json">(.*?)</script>', html, re.S)
    embedded = match.group(1)
    assert json.loads(embedded) == result
    assert "<" not in embedded and "&" not in embedded
    assert "\u2028" not in embedded and "\u2029" not in embedded
    assert "\\u2028" in embedded and "\\u2029" in embedded
    assert ".innerHTML" not in html and "insertAdjacentHTML" not in html


def test_nonfinite_data_rejected():
    with pytest.raises(ValueError):
        _script_json({"score": float("nan")})


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_browser_constraints_preserve_phase_and_detect_conflicts():
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    test = r"""
const assert = require('node:assert/strict');
const core = require(process.argv[1]);
const scenarios = [
  { id: 'S1', orientations: {A: 0, B: 0, C: 0}, outcome_key: 'one' },
  { id: 'S2', orientations: {A: 0, B: 0, C: 1}, outcome_key: 'two' },
  { id: 'S3', orientations: {A: 0, B: 1, C: 0}, outcome_key: 'two' },
  { id: 'S4', orientations: {A: 0, B: 1, C: 1}, outcome_key: 'three' }
];
const snapshot = JSON.stringify(scenarios);
assert.equal(core.filterScenarios(scenarios, []).length, 4);
assert.deepEqual(core.filterScenarios(scenarios, [{a:'A',b:'B',parity:1}]).map(x=>x.id), ['S3','S4']);
assert.deepEqual(core.filterScenarios(scenarios, [{a:'B',b:'A',parity:1}]).map(x=>x.id), ['S3','S4']);
assert.equal(core.filterScenarios(scenarios, [{a:'A',b:'B',parity:0},{a:'B',b:'C',parity:0},{a:'A',b:'C',parity:1}]).length, 0);
assert.equal(core.filterScenarios(scenarios, [{a:'A',b:'missing',parity:0}]).length, 0);
assert.deepEqual(core.groupScenarios(scenarios).map(x=>x.length), [1,2,1]);
assert.equal(core.filterScenarios([{id:'fixed',orientations:{},outcome_key:'fixed'}], []).length, 1);
assert.equal(JSON.stringify(scenarios), snapshot);
"""
    subprocess.run(["node", "-e", test, str(javascript)], check=True, capture_output=True, text=True)
