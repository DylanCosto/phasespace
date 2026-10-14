"""Biological fixtures exercise joint edits, phase constraints and edge cases."""

from dataclasses import replace
import json

import pytest

from phasespace.analysis import analyze
from phasespace.models import Evidence, Locus, PhaseSpaceError, Transcript, Variant


def variant(
    identifier,
    start,
    ref,
    alt,
    *,
    gt=(0, 1),
    phased=False,
    ps=None,
    alleles=None,
):
    called = alleles if alleles is not None else tuple(ref if g == 0 else alt for g in gt)
    return Variant(
        identifier, start + 100, ref, called, gt, phased, ps,
        start, start + len(ref), ref, called, 60, 20,
    )


def locus(reference, variants=(), *, strand="+"):
    transcript = Transcript("tx1", "SYNTHETIC", "chr1", strand, ((100, 100 + len(reference)),), reference)
    return Locus(transcript, "sample1", list(variants), [Evidence(100, 100 + len(reference), "unknown")])


def changed_haplotype(result):
    return next(h for h in result["scenarios"][0]["haplotypes"] if h["variant_ids"])


def test_cis_trans_stops_keep_both_biological_possibilities():
    data = locus("ATGCAACAAGGGTAA", [variant("A", 3, "C", "T"), variant("B", 6, "C", "T")])
    result = analyze(data)
    assert result["scenario_count"] == result["outcome_count"] == 2
    assert result["reference_protein"] == "MQQG"
    assert result["reference_first_stop"] == 5
    cis, trans = result["scenarios"]
    assert cis["unaffected_copies"] == 1
    assert [h["first_stop"] for h in cis["haplotypes"]] == [5, 2]
    assert cis["haplotypes"][1]["variant_ids"] == ["A", "B"]
    assert trans["unaffected_copies"] == 0
    assert [h["first_stop"] for h in trans["haplotypes"]] == [3, 2]
    assert all(h["consequences"] == ["stop_gained"] for h in trans["haplotypes"])
    assert result["informative_pairs"][0]["score"] == 1
    assert json.loads(json.dumps(result))["schema_version"] == 1


@pytest.mark.parametrize("ps", ["100", None])
def test_phased_gt_orientation_is_preserved_in_explicit_and_implicit_blocks(ps):
    data = locus("ATGCAACAAGGGTAA", [
        variant("A", 3, "C", "T", gt=(1, 0), phased=True, ps=ps),
        variant("B", 6, "C", "T", gt=(0, 1), phased=True, ps=ps),
    ])
    result = analyze(data)
    assert len(result["blocks"]) == 1
    assert result["scenario_count"] == 1
    scenario = result["scenarios"][0]
    assert [h["variant_ids"] for h in scenario["haplotypes"]] == [["A"], ["B"]]
    assert scenario["unaffected_copies"] == 0
    assert any("implicit phase set" in w for w in result["warnings"]) == (ps is None)


def test_distinct_phase_sets_are_independent_and_unphased_ps_is_not_used():
    result = analyze(locus("ATGCAACAACAATAA", [
        variant("A", 3, "C", "T", phased=True, ps="a"),
        variant("B", 6, "C", "T", phased=True, ps="b"),
        variant("C", 9, "C", "T", phased=False, ps="a"),
    ]))
    assert len(result["blocks"]) == 3
    assert result["scenario_count"] == 4
    assert all(s["orientations"]["block-1"] == 0 for s in result["scenarios"])


def test_homozygote_is_fixed_on_both_copies_and_needs_no_block():
    result = analyze(locus("ATGCAATAA", [variant("A", 3, "C", "T", gt=(1, 1))]), max_blocks=0)
    assert result["blocks"] == []
    assert result["variants"][0]["block"] is None
    assert result["scenario_count"] == 1
    assert result["scenarios"][0]["orientations"] == {}
    assert result["scenarios"][0]["unaffected_copies"] == 0
    assert [h["protein"] for h in result["scenarios"][0]["haplotypes"]] == ["M", "M"]


def test_no_variants_yields_one_reference_scenario():
    result = analyze(locus("ATGCAATAA"), max_blocks=0)
    assert result["scenario_count"] == result["outcome_count"] == 1
    assert result["scenarios"][0]["unaffected_copies"] == 2
    assert result["informative_pairs"] == []
    assert all(h["consequences"] == ["reference"] for h in result["scenarios"][0]["haplotypes"])


def test_multiallelic_one_two_uses_both_called_alleles():
    result = analyze(locus("ATGCAATAA", [
        variant("A", 3, "C", "T", gt=(1, 2), alleles=("T", "A")),
    ]))
    copies = result["scenarios"][0]["haplotypes"]
    assert [h["protein"] for h in copies] == ["M", "MK"]
    assert [h["variant_ids"] for h in copies] == [["A"], ["A"]]
    assert result["scenarios"][0]["unaffected_copies"] == 0


def test_joint_codon_edits_remove_single_variant_stop_prediction():
    result = analyze(locus("ATGCAATAA", [
        variant("A", 3, "C", "T", phased=True, ps="one"),
        variant("B", 4, "A", "C", phased=True, ps="one"),
    ]))
    copy = changed_haplotype(result)
    assert copy["cds"] == "ATGTCATAA"
    assert copy["protein"] == "MS"
    assert copy["first_stop"] == 3
    assert copy["consequences"] == ["missense"]


def test_compensating_indels_preserve_intervening_protein_changes():
    edits = [
        variant("ins", 3, "A", "AT", phased=True, ps="one"),
        variant("del", 9, "GG", "G", phased=True, ps="one"),
    ]
    result = analyze(locus("ATGAAACCCGGGTAA", edits))
    copy = changed_haplotype(result)
    assert copy["cds"] == "ATGATAACCCGGTAA"
    assert copy["protein"] == "MITR"
    assert result["reference_protein"] == "MKPG"
    assert "frame_restored" in copy["consequences"]
    assert "frameshift" not in copy["consequences"]
    assert "stop_gained" not in copy["consequences"]
    assert result["scenarios"][0]["unaffected_copies"] == 1
    reversed_result = analyze(locus("ATGAAACCCGGGTAA", list(reversed(edits))))
    assert result["scenarios"] == reversed_result["scenarios"]


@pytest.mark.parametrize("ref,alt,expected", [("GAAA", "G", "inframe_deletion"), ("G", "GAAA", "inframe_insertion")])
def test_moving_terminal_stop_with_inframe_edit_does_not_gain_or_lose_it(ref, alt, expected):
    result = analyze(locus("ATGAAACCCTAA", [variant("indel", 2, ref, alt)]))
    copy = changed_haplotype(result)
    assert expected in copy["consequences"]
    assert "stop_gained" not in copy["consequences"]
    assert "stop_lost" not in copy["consequences"]


def test_single_base_indel_reports_frameshift_and_translates_complete_codons():
    copy = changed_haplotype(analyze(locus("ATGAAACCCTAA", [variant("ins", 3, "A", "AT")])))
    assert copy["cds"] == "ATGATAACCCTAA"
    assert copy["protein"] == "MITL"
    assert copy["protein_length"] == 4
    assert copy["first_stop"] is None
    assert "frameshift" in copy["consequences"]
    assert "stop_lost" in copy["consequences"]


def test_frameshift_with_premature_stop_does_not_also_report_stop_lost():
    copy = changed_haplotype(analyze(locus("ATGCTAAGCCCCTAA", [variant("del", 2, "GC", "G")])))
    assert copy["cds"] == "ATGTAAGCCCCTAA"
    assert copy["first_stop"] == 2
    assert copy["consequences"] == ["frameshift", "stop_gained"]


def test_negative_strand_uses_preprojected_cds_alleles_without_reversing_again():
    edit = replace(variant("A", 3, "C", "T"), pos=104, ref="G", alleles=("G", "A"))
    result = analyze(locus("ATGCAATAA", [edit], strand="-"))
    assert changed_haplotype(result)["cds"] == "ATGTAATAA"
    assert changed_haplotype(result)["protein"] == "M"


def test_synonymous_does_not_count_as_reference_identical():
    result = analyze(locus("ATGCAATAA", [variant("syn", 5, "A", "G")]))
    assert changed_haplotype(result)["consequences"] == ["synonymous"]
    assert result["scenarios"][0]["unaffected_copies"] == 1


def test_start_lost_and_stop_lost_are_explicit():
    result = analyze(locus("ATGCAATAA", [
        variant("start", 0, "A", "G", gt=(1, 1)),
        variant("stop", 6, "T", "C", gt=(1, 1)),
    ]))
    copy = result["scenarios"][0]["haplotypes"][0]
    assert copy["protein"] == "VQQ"
    assert {"start_lost", "stop_lost", "missense"} <= set(copy["consequences"])


def test_downstream_missense_is_not_translated_after_a_gained_stop():
    copy = changed_haplotype(analyze(locus("ATGCAACAATAA", [
        variant("stop", 3, "C", "T", phased=True, ps="one"),
        variant("downstream", 6, "C", "A", phased=True, ps="one"),
    ])))
    assert copy["protein"] == "M"
    assert copy["consequences"] == ["stop_gained"]


def test_ambiguous_codon_is_not_called_missense_or_definite_stop():
    result = analyze(locus("ATGCAATAA", [variant("uncertain", 5, "A", "N")]))
    copy = changed_haplotype(result)
    assert copy["protein"] == "MX"
    assert copy["consequences"] == ["uncertain_codon"]
    assert any("N-containing codons" in w for w in result["warnings"])


def test_resolving_uncertain_reference_codon_to_stop_does_not_assert_stop_gain():
    copy = changed_haplotype(analyze(locus("ATGTNATAA", [variant("resolved", 4, "N", "A")])))
    assert copy["protein"] == "M"
    assert copy["first_stop"] == 2
    assert copy["consequences"] == ["uncertain_codon"]


def test_outcomes_use_full_cds_not_only_translated_protein():
    result = analyze(locus("ATGCAACAATAA", [variant("synA", 5, "A", "G"), variant("synB", 8, "A", "G")]))
    assert result["scenario_count"] == result["outcome_count"] == 2
    assert result["protein_outcome_count"] == 1
    assert result["informative_pairs"][0]["score"] == 0
    assert all(h["protein"] == "MQQ" for s in result["scenarios"] for h in s["haplotypes"])


def test_an_existing_later_stop_is_not_called_gained_after_stop_loss():
    copy = changed_haplotype(analyze(locus("ATGTAACAATAA", [variant("stop", 3, "T", "C")])))
    assert copy["protein"] == "MQQ"
    assert copy["first_stop"] == 4
    assert copy["consequences"] == ["stop_lost"]


@pytest.mark.parametrize("position,ref", [(0, "A"), (6, "T")])
def test_unknown_start_or_stop_base_does_not_assert_loss(position, ref):
    copy = changed_haplotype(analyze(locus("ATGCAATAA", [variant("unknown", position, ref, "N")])))
    assert "uncertain_codon" in copy["consequences"]
    assert "start_lost" not in copy["consequences"]
    assert "stop_lost" not in copy["consequences"]


def test_global_copy_swap_does_not_change_outcome_key():
    first = variant("A", 3, "C", "T", phased=True, ps="p")
    second = variant("B", 6, "C", "T", phased=True, ps="p")
    normal = analyze(locus("ATGCAACAATAA", [first, second]))
    swapped = analyze(locus("ATGCAACAATAA", [
        replace(v, gt=(1, 0), alleles=tuple(reversed(v.alleles)), cds_alleles=tuple(reversed(v.cds_alleles)))
        for v in [first, second]
    ]))
    assert normal["scenarios"][0]["outcome_key"] == swapped["scenarios"][0]["outcome_key"]


def test_equivalent_outcomes_keep_arrangements_and_have_no_informative_pair():
    # The two ALT indexes denote identical supplied coding sequences. Genotype
    # arrangements are retained, but their biological sequence outcome is one.
    result = analyze(locus("ATGCAACAATAA", [
        variant("same_sequence", 3, "C", "T", gt=(1, 2), alleles=("T", "T")),
        variant("B", 6, "C", "T"),
    ]))
    assert result["scenario_count"] == 2
    assert result["outcome_count"] == 1
    assert result["informative_pairs"][0]["score"] == 0


def test_excess_blocks_raise_without_partial_results():
    data = locus("ATGCAACAACAATAA", [variant(str(i), i, "C", "T") for i in [3, 6, 9]])
    with pytest.raises(PhaseSpaceError, match="3 independent phase blocks.*max_blocks=2"):
        analyze(data, max_blocks=2)


@pytest.mark.parametrize("limit", [-1, 13, 1.5, True])
def test_invalid_limits_are_rejected(limit):
    with pytest.raises(PhaseSpaceError, match="max_blocks"):
        analyze(locus("ATGCAATAA"), max_blocks=limit)


def test_overlaps_are_rejected_even_if_an_arrangement_could_avoid_them():
    data = locus("ATGAAATAA", [variant("A", 3, "AA", "A"), variant("B", 4, "A", "C")])
    with pytest.raises(PhaseSpaceError, match="Overlapping coding edits"):
        analyze(data)


@pytest.mark.parametrize("reference", ["ATGA", "", "ATGRRRTAA"])
def test_invalid_reference_is_rejected(reference):
    with pytest.raises(PhaseSpaceError, match="reference CDS"):
        analyze(locus(reference))


@pytest.mark.parametrize("change,match", [
    ({"cds_ref": "G"}, "reference allele does not match"),
    ({"cds_start": -1}, "outside the reference CDS"),
    ({"gt": (None, 1)}, "nonnegative integers"),
    ({"cds_alleles": ("T", "T")}, "GT=0"),
    ({"gt": (1, 1)}, "identical genotype indexes"),
    ({"cds_alleles": ("C", "<DEL>")}, "unsupported symbolic"),
])
def test_inconsistent_alleles_are_rejected(change, match):
    edit = replace(variant("A", 3, "C", "T"), **change)
    with pytest.raises(PhaseSpaceError, match=match):
        analyze(locus("ATGCAATAA", [edit]))


def test_duplicate_identifiers_are_rejected():
    with pytest.raises(PhaseSpaceError, match="identifiers"):
        analyze(locus("ATGCAACAATAA", [variant("A", 3, "C", "T"), variant("A", 6, "C", "T")]))


def test_missing_evidence_warning_is_conditional_and_not_duplicated():
    data = locus("ATGCAATAA")
    data.warnings = ["Some CDS bases lack usable genotype evidence; additional variants may be unobserved."]
    assert analyze(data)["warnings"] == data.warnings
    data.evidence = [Evidence(100, 109, "reference")]
    data.warnings = []
    assert analyze(data)["warnings"] == []


def test_annotation_without_terminal_stop_gets_concise_notice():
    result = analyze(locus("ATGCAA"))
    assert result["reference_first_stop"] is None
    assert result["reference_protein"] == "MQ"
    assert any("no confirmed terminal stop" in w for w in result["warnings"])
