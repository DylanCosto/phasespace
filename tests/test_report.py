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


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_protein_scale_stays_comparable_across_arrangements():
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    test = r"""
const assert = require('node:assert/strict');
const { sharedProteinScale } = require(process.argv[1]);
const scenarios = [
  { id: 'short', haplotypes: [{protein_length: 7}, {protein_length: 20}] },
  { id: 'long', haplotypes: [{protein_length: 31}, {protein_length: 45}] },
  { id: 'empty-translation', haplotypes: [{protein_length: 0}, {protein_length: 0}] }
];
const snapshot = JSON.stringify(scenarios);
// A longer product in another arrangement must set the scale for every copy.
assert.equal(sharedProteinScale(scenarios, 31), 45);
assert.equal(sharedProteinScale([...scenarios].reverse(), 31), 45);
// Reference length remains the lower bound when all predictions are shorter.
assert.equal(sharedProteinScale(scenarios, 60), 60);
assert.equal(sharedProteinScale([], 31), 31);
// Zero-length and absent translations still need a nonzero plotting denominator.
assert.equal(sharedProteinScale([], 0), 1);
assert.equal(sharedProteinScale([], undefined), 1);
assert.equal(sharedProteinScale([scenarios[2]], 0), 1);
assert.equal(JSON.stringify(scenarios), snapshot);
"""
    subprocess.run(["node", "-e", test, str(javascript)], check=True, capture_output=True, text=True)


def _embedded(html):
    match = re.search(r'<script id="phasespace-data" type="application/json">(.*?)</script>', html, re.S)
    return json.loads(match.group(1))


def _joint_edit_result():
    from phasespace.analysis import analyze
    from phasespace.models import Evidence, Locus, Transcript, Variant

    reference = "ATGCAACAAGGGAAACCCNTTGGGTAA"

    def call(identifier, start, ref, alleles, gt, phased=False, ps=None):
        return Variant(identifier, start + 100, ref, alleles, gt, phased, ps, start, start + len(ref), ref, alleles, 60, 20)

    variants = [
        call("stop", 3, "C", ("C", "T"), (0, 1)),
        call("ins", 8, "A", ("AT", "A"), (1, 0), True, "7"),
        call("del", 13, "AA", ("AA", "A"), (0, 1), True, "7"),
        call("multi", 16, "C", ("G", "T"), (1, 2)),
        call("hom", 22, "G", ("A", "A"), (1, 1)),
    ]
    transcript = Transcript("tx1", "GENE", "chr1", "+", ((100, 100 + len(reference)),), reference)
    return analyze(Locus(transcript, "sample1", variants, [Evidence(100, 127, "unknown")]))


def test_report_embeds_edits_instead_of_full_sequences(tmp_path):
    result = _joint_edit_result()
    output = tmp_path / "report.html"
    write_report(result, output)
    embedded = _embedded(output.read_text())
    assert embedded["report_encoding"] == "edits"
    copies = [copy for scenario in embedded["scenarios"] for copy in scenario["haplotypes"]]
    assert copies and all("cds" not in copy and "protein" not in copy and "alleles" in copy for copy in copies)
    # Nothing else in the analysis changes.
    full = json.loads(json.dumps(result))
    for key in full:
        if key != "scenarios":
            assert embedded[key] == full[key]


def test_compact_report_size_does_not_scale_with_sequence_per_arrangement(tmp_path):
    from phasespace.analysis import analyze
    from phasespace.models import Evidence, Locus, Transcript, Variant

    reference = "ATG" + "GCT" * 1998 + "TAA"
    variants = [Variant(str(i), 10 + 300 * i, "C", ("C", "A"), (0, 1), False, None, 10 + 300 * i, 11 + 300 * i, "C", ("C", "A"), 60, 30)
                for i in range(8)]
    transcript = Transcript("tx", "G", "c", "+", ((0, len(reference)),), reference)
    result = analyze(Locus(transcript, "s", variants, [Evidence(0, len(reference), "reference")]))
    output = tmp_path / "report.html"
    write_report(result, output)
    assert len(json.dumps(result)) > 2_000_000
    assert output.stat().st_size < 250_000


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_browser_rebuilds_exact_analysis_from_compact_report(tmp_path):
    result = _joint_edit_result()
    output = tmp_path / "report.html"
    write_report(result, output)
    compact = tmp_path / "compact.json"
    compact.write_text(json.dumps(_embedded(output.read_text())))
    expected = tmp_path / "expected.json"
    expected.write_text(json.dumps(result))
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    test = r"""
const assert = require('node:assert/strict');
const fs = require('node:fs');
const core = require(process.argv[1]);
const compact = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const expected = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const snapshot = JSON.stringify(compact);
assert.deepEqual(core.expandAnalysis(compact), expected);
assert.equal(JSON.stringify(compact), snapshot);
// Translation matches the engine: N codons become X and stop ends the protein.
assert.equal(core.translate('ATGNAATAAGGG'), 'MX');
assert.equal(core.translate('ATGAA'), 'M');
"""
    subprocess.run(["node", "-e", test, str(javascript), str(compact), str(expected)], check=True, capture_output=True, text=True)


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_browser_verdict_flip_and_change_labels():
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    test = r"""
const assert = require('node:assert/strict');
const core = require(process.argv[1]);
const ref = {consequences: ['reference']};
const alt = {consequences: ['stop_gained']};
const cis = {id: 'cis', orientations: {A: 0, B: 0, C: 0}, haplotypes: [ref, alt]};
const trans = {id: 'trans', orientations: {A: 0, B: 1, C: 0}, haplotypes: [alt, alt]};
assert.equal(core.phaseVerdict([cis], 'X').tone, 'good');
assert.equal(core.phaseVerdict([trans], 'X').tone, 'none');
assert.deepEqual(core.phaseVerdict([cis, trans], 'X'), {tone: 'mixed', total: 2, keep: 1, min: 0, max: 1, cdsChanged: true});
// Synonymous edits change the CDS but keep the reference protein.
const syn = {consequences: ['synonymous']};
const both = core.phaseVerdict([{haplotypes: [syn, syn]}, {haplotypes: [ref, syn]}], 'X');
assert.equal(both.tone, 'good'); assert.equal(both.min, 2); assert.equal(both.cdsChanged, true);
assert.equal(core.phaseVerdict([{haplotypes: [ref, ref]}], 'X').cdsChanged, false);
assert.equal(core.phaseVerdict([], 'X').tone, 'empty');
// Older results without consequence labels compare full sequences.
assert.equal(core.phaseVerdict([{haplotypes: [{cds: 'X'}, {cds: 'Y'}]}], 'X').tone, 'good');

const all = [0, 1].flatMap(b => [0, 1].map(c => ({id: `${b}${c}`, orientations: {A: 0, B: b, C: c}})));
assert.equal(core.flipTarget(all, ['A', 'B', 'C'], all[0], 'B').id, '10');
assert.equal(core.flipTarget(all, ['A', 'B', 'C'], all[0], 'C').id, '01');
// The first block is fixed: flipping it is the same as flipping every other block.
assert.equal(core.flipTarget(all, ['A', 'B', 'C'], all[0], 'A').id, '11');
assert.equal(core.flipTarget([all[0]], ['A'], all[0], 'A'), null);

const v = (start, ref) => ({cds_start: start, cds_end: start + ref.length, cds_ref: ref});
assert.equal(core.changeLabel(v(21, 'C'), 'T'), 'c.22C>T');
assert.equal(core.changeLabel(v(30, 'CCA'), 'C'), 'c.32_33del');
assert.equal(core.changeLabel(v(30, 'CA'), 'C'), 'c.32del');
assert.equal(core.changeLabel(v(30, 'GTC'), 'C'), 'c.31_32del');
assert.equal(core.changeLabel(v(8, 'A'), 'AT'), 'c.9_10insT');
assert.equal(core.changeLabel(v(8, 'A'), 'TA'), 'c.8_9insT');
assert.equal(core.changeLabel(v(4, 'AC'), 'GT'), 'c.5_6delinsGT');
"""
    subprocess.run(["node", "-e", test, str(javascript)], check=True, capture_output=True, text=True)


@pytest.mark.skipif(shutil.which("node") is None, reason="Node is needed for browser logic test")
def test_pin_labels_never_overlap_or_cross_stems():
    javascript = Path(__file__).parents[1] / "src" / "phasespace" / "assets" / "report.js"
    test = r"""
const assert = require('node:assert/strict');
const { placeLabels } = require(process.argv[1]);
// A dense real-world cluster: 19 labels within 20% of the track, plus spread-out ones.
const items = [];
for (let i = 0; i < 19; i++) items.push({center: 200 + i * 9, width: 70});
items.push({center: 50, width: 70}, {center: 900, width: 70});
items.sort((a, b) => a.center - b.center);
const spots = placeLabels(items, 1000, 3);
const placed = spots.map((s, i) => ({...s, center: items[i].center, right: s.left + items[i].width})).filter(s => s.row >= 0);
assert.ok(placed.every(s => s.row < 3));
for (const a of placed) for (const b of placed) {
  if (a === b) continue;
  if (a.row === b.row) assert.ok(a.right + 6 <= b.left || b.right + 6 <= a.left, 'labels overlap');
  if (a.row > b.row) assert.ok(a.center <= b.left - 3 || a.center >= b.right + 3, 'stem crosses label');
}
// Isolated labels always get the first row; the dense cluster cannot all fit.
assert.equal(spots[0].row, 0); assert.equal(spots[spots.length - 1].row, 0);
assert.ok(spots.some(s => s.row === -1));
// Labels are clamped inside the track.
assert.ok(spots.every((s, i) => s.left >= 0 && s.left + items[i].width <= 1000));
"""
    subprocess.run(["node", "-e", test, str(javascript)], check=True, capture_output=True, text=True)
