"""Biological input guarantees for sequence projection and genotype evidence."""
from pathlib import Path

import pysam
import pytest

from phasespace.inputs import load_locus
from phasespace.models import PhaseSpaceError


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=90>
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum block depth">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""


@pytest.fixture
def files(tmp_path):
    # Coding intervals [12,18), [30,39): ATGAAA + CCCGGGTAA.
    sequence = list("A" * 90)
    sequence[12:18] = "ATGAAA"
    sequence[30:39] = "CCCGGGTAA"
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(reference))
    annotation = tmp_path / "genes.gtf"
    annotation.write_text(
        'chr1\ttest\tCDS\t13\t18\t.\t+\t0\tgene_id "GENE"; transcript_id "TX"; gene_name "Gene";\n'
        'chr1\ttest\tCDS\t31\t39\t.\t+\t0\tgene_id "GENE"; transcript_id "TX"; gene_name "Gene";\n'
    )
    vcf = tmp_path / "calls.vcf"

    def write(records, *, indexed=False, header=HEADER):
        vcf.write_text(header + "\n".join(records) + ("\n" if records else ""))
        if indexed:
            compressed = tmp_path / "calls.vcf.gz"
            pysam.tabix_compress(str(vcf), str(compressed), force=True)
            pysam.tabix_index(str(compressed), preset="vcf", force=True)
            return compressed
        return vcf

    return reference, annotation, write


def load(files, records, **kwargs):
    reference, annotation, write = files
    return load_locus(write(records), reference, annotation, "TX", **kwargs)


def row(pos, ref, alt, gt="0/1", *, info=".", filt="PASS", fmt="GT:GQ:DP", values=None):
    values = values or f"{gt}:50:12"
    return f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\t{fmt}\t{values}"


def test_positive_strand_multiallelic_gt_order_and_phase_set(files):
    locus = load(files, [row(16, "A", "T,G", fmt="GT:GQ:DP:PS", values="2|1:50:12:13")])
    assert locus.transcript.reference_cds == "ATGAAACCCGGGTAA"
    assert locus.transcript.gene == "Gene"
    variant = locus.variants[0]
    assert variant.gt == (2, 1)
    assert variant.alleles == ("G", "T")
    assert variant.cds_alleles == ("G", "T")
    assert (variant.cds_start, variant.cds_end, variant.cds_ref) == (3, 4, "A")
    assert variant.phased and variant.phase_set == "13"
    assert any(e.status == "unknown" for e in locus.evidence)
    assert any("streamed" in warning for warning in locus.warnings)


def test_negative_strand_reverse_complements_alleles_and_cds_offsets(files):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text().replace("\t+\t", "\t-\t"))
    locus = load_locus(write([row(32, "C", "CT", fmt="GT:GQ:DP:PS", values="1|0:50:12:7")]), reference, annotation, "TX")
    assert locus.transcript.reference_cds == "TTACCCGGGTTTCAT"
    variant = locus.variants[0]
    assert variant.alleles == ("CT", "C")
    assert variant.cds_alleles == ("AG", "G")
    assert (variant.cds_start, variant.cds_end, variant.cds_ref) == (7, 8, "G")
    assert variant.gt == (1, 0) and variant.phase_set == "7"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_internal_exon_phase_bases_are_retained(files, strand):
    reference, annotation, write = files
    # Length 4 then 8 on +; length 8 then 4 on -. Both join to 12 bases.
    first_phase, second_phase = (0, 2) if strand == "+" else (1, 0)
    annotation.write_text(
        f'chr1\ttest\tCDS\t13\t16\t.\t{strand}\t{first_phase}\ttranscript_id "TX";\n'
        f'chr1\ttest\tCDS\t31\t38\t.\t{strand}\t{second_phase}\ttranscript_id "TX";\n'
    )
    locus = load_locus(write([]), reference, annotation, "TX")
    assert len(locus.transcript.reference_cds) == 12


def test_gff3_gene_resolution_and_multiple_parents(files):
    reference, annotation, write = files
    annotation.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene1;Name=Example%20Gene\n"
        "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=TX;Parent=gene1\n"
        "chr1\ttest\tCDS\t13\t18\t.\t+\t0\tParent=OTHER,TX\n"
        "chr1\ttest\tCDS\t31\t39\t.\t+\t0\tParent=TX\n"
    )
    locus = load_locus(write([]), reference, annotation, "TX")
    assert locus.transcript.gene == "Example Gene"


@pytest.mark.parametrize("indexed", [False, True])
def test_reference_block_beginning_before_query_covers_only_cds(files, indexed):
    reference, annotation, write = files
    path = write([row(1, "A", "<NON_REF>", info="END=50", fmt="GT:GQ:DP:MIN_DP", values="0/0:40:20:9")], indexed=indexed)
    locus = load_locus(path, reference, annotation, "TX")
    assert locus.variants == []
    assert [(e.start, e.end, e.status) for e in locus.evidence] == [(12, 18, "reference"), (30, 39, "reference")]
    if indexed:
        assert not any("streamed" in warning for warning in locus.warnings)


def test_block_min_depth_not_anchor_depth_drives_evidence(files):
    locus = load(files, [row(1, "A", "<NON_REF>", info="END=50", fmt="GT:GQ:DP:MIN_DP", values="0/0:40:20:2")])
    assert all(item.status == "unknown" for item in locus.evidence)
    assert any("depth below 5" in warning for warning in locus.warnings)


def test_single_base_reference_block_can_use_min_dp(files):
    locus = load(files, [row(16, "A", "<NON_REF>", info="END=16", fmt="GT:GQ:MIN_DP", values="0/0:40:9")])
    assert [(item.start, item.end) for item in locus.evidence if item.status == "reference"] == [(15, 16)]


@pytest.mark.parametrize("fmt,values,filt,reason", [
    ("GT:DP", "0/0:20", "PASS", "missing GQ"),
    ("GT:GQ", "0/0:50", "PASS", "missing MIN_DP"),
    ("GT:GQ:DP", "0/0:50:20", "LowQual", "FILTER=LowQual"),
    ("GT:GQ:DP:FT", "0/0:50:20:LowQual", "PASS", "FT=LowQual"),
    ("GT:GQ:DP", "0/.:50:20", "PASS", "partial genotype"),
    ("GT:GQ:DP", "./.:50:20", "PASS", "partial genotype"),
])
def test_unreliable_reference_block_does_not_assert_reference(files, fmt, values, filt, reason):
    locus = load(files, [row(1, "A", "<NON_REF>", info="END=50", fmt=fmt, values=values, filt=filt)])
    assert all(item.status == "unknown" for item in locus.evidence)
    assert any(reason in warning for warning in locus.warnings)


def test_missing_variant_quality_excludes_call_and_preserves_unknown(files):
    locus = load(files, [row(16, "A", "T", fmt="GT:DP", values="0/1:12")])
    assert locus.variants == []
    assert all(item.status == "unknown" for item in locus.evidence)
    assert any("missing GQ" in warning for warning in locus.warnings)


def test_partial_alt_call_does_not_become_a_phased_variant(files):
    locus = load(files, [row(16, "A", "T", fmt="GT:GQ:DP", values="1/.:50:12")])
    assert locus.variants == []
    assert all(item.status == "unknown" for item in locus.evidence)


def test_plain_vcf_absence_stays_unknown(files):
    locus = load(files, [])
    assert [(e.start, e.end, e.status) for e in locus.evidence] == [(12, 18, "unknown"), (30, 39, "unknown")]


def test_homozygous_reference_record_asserts_only_its_ref_span(files):
    locus = load(files, [row(16, "AAA", "A", gt="0/0")])
    assert [(e.start, e.end, e.status) for e in locus.evidence] == [(12, 15, "unknown"), (15, 18, "reference"), (30, 39, "unknown")]


def test_full_ref_validation_checks_bases_outside_cds_before_boundary_error(files):
    # The REF anchor is outside CDS, and its second base is wrong.
    with pytest.raises(PhaseSpaceError, match="REF mismatch"):
        load(files, [row(12, "AC", "A")])


def test_full_ref_validation_even_when_call_low_quality(files):
    with pytest.raises(PhaseSpaceError, match="REF mismatch"):
        load(files, [row(16, "C", "T", values="0/1:0:0")])


def test_called_variant_crossing_cds_boundary_is_explicit_error(files):
    with pytest.raises(PhaseSpaceError, match="CDS boundary"):
        load(files, [row(12, "AA", "A")])


def test_insertion_at_cds_edge_is_explicit_error(files):
    with pytest.raises(PhaseSpaceError, match="CDS boundary"):
        load(files, [row(18, "A", "AT")])


@pytest.mark.parametrize("alt,info", [("<DEL>", "END=40"), ("*", "."), ("A]chr1:50]", "."), ("<NON_REF>", ".")])
def test_unsupported_called_alleles_rejected_even_after_supported_call(files, alt, info):
    with pytest.raises(PhaseSpaceError, match="Unsupported called allele"):
        load(files, [row(16, "A", "T"), row(32, "C", alt, info=info)])


def test_unused_non_ref_allele_is_allowed(files):
    locus = load(files, [row(16, "A", "T,<NON_REF>", gt="0/1")])
    assert locus.variants[0].alleles == ("A", "T")


@pytest.mark.parametrize("genotype", ["1", "0/1/1"])
def test_non_diploid_genotype_explicitly_rejected(files, genotype):
    with pytest.raises(PhaseSpaceError, match="Non-diploid"):
        load(files, [row(16, "A", "T"), row(32, "C", "T", gt=genotype)])


def test_overlapping_variant_edits_explicitly_rejected(files):
    with pytest.raises(PhaseSpaceError, match="Overlapping"):
        load(files, [row(16, "AA", "A"), row(17, "A", "T")])


def test_unphased_call_does_not_inherit_ps(files):
    locus = load(files, [row(16, "A", "T", fmt="GT:GQ:DP:PS", values="0/1:50:12:13")])
    assert locus.variants[0].phase_set is None
    assert not locus.variants[0].phased


def test_phased_call_without_ps_preserves_phased_flag(files):
    locus = load(files, [row(16, "A", "T", gt="1|0")])
    assert locus.variants[0].phased
    assert locus.variants[0].gt == (1, 0)
    assert locus.variants[0].phase_set is None


def test_inconsistent_cds_phase_rejected(files):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text().replace("31\t39\t.\t+\t0", "31\t39\t.\t+\t1"))
    with pytest.raises(PhaseSpaceError, match="Inconsistent CDS phase"):
        load_locus(write([]), reference, annotation, "TX")


def test_partial_cds_rejected(files):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text().replace("13\t18\t.\t+\t0", "13\t18\t.\t+\t1"))
    with pytest.raises(PhaseSpaceError, match="Partial CDS"):
        load_locus(write([]), reference, annotation, "TX")


def test_reference_index_required_without_writing_user_reference(files):
    reference, annotation, write = files
    Path(str(reference) + ".fai").unlink()
    with pytest.raises(PhaseSpaceError, match="existing .fai"):
        load_locus(write([]), reference, annotation, "TX")
    assert not Path(str(reference) + ".fai").exists()


def test_sample_must_be_selected_for_cohort(files):
    reference, annotation, write = files
    header = HEADER.replace("\tSAMPLE\n", "\tSAMPLE\tSECOND\n")
    path = write([], header=header)
    with pytest.raises(PhaseSpaceError, match="Choose a sample"):
        load_locus(path, reference, annotation, "TX")
    assert load_locus(path, reference, annotation, "TX", sample="SECOND").sample == "SECOND"


def test_overlapping_passing_reference_and_variant_calls_rejected(files):
    with pytest.raises(PhaseSpaceError, match="Conflicting passing reference and variant"):
        load(files, [row(1, "A", "<NON_REF>", info="END=50", fmt="GT:GQ:DP:MIN_DP", values="0/0:50:12:10"), row(16, "A", "T")])


def ref_block(pos, end, base="A"):
    return row(pos, base, "<NON_REF>", info=f"END={end}", fmt="GT:GQ:DP:MIN_DP", values="0/0:50:12:10")


def test_gvcf_reference_block_through_deleted_bases_is_not_a_conflict(files):
    # Some gVCF writers resume reference blocks inside a heterozygous deletion.
    locus = load(files, [ref_block(13, 30), row(31, "CCC", "C"), ref_block(32, 39, "C")])
    assert [v.cds_alleles for v in locus.variants] == [("CCC", "C")]
    spans = [(e.start, e.end, e.status) for e in locus.evidence]
    assert spans == [(12, 18, "reference"), (30, 33, "variant"), (33, 39, "reference")]


def test_reference_block_over_deletion_anchor_still_conflicts(files):
    with pytest.raises(PhaseSpaceError, match="Conflicting passing reference and variant"):
        load(files, [ref_block(13, 31), row(31, "CCC", "C"), ref_block(32, 39, "C")])


def test_spanning_deletion_allele_covered_by_upstream_deletion_is_accepted(files):
    locus = load(files, [row(31, "CCC", "C"), row(32, "C", "*"), row(33, "C", "*,<NON_REF>")])
    assert [v.cds_alleles for v in locus.variants] == [("CCC", "C")]
    variant_bases = {base for e in locus.evidence if e.status == "variant" for base in range(e.start, e.end)}
    assert variant_bases == {30, 31, 32}


def test_low_quality_spanning_deletion_record_stays_unknown(files):
    locus = load(files, [row(31, "CCC", "C"), row(32, "C", "*", values="0/1:5:12")])
    assert any(e.status == "unknown" and e.start == 31 for e in locus.evidence)


def test_spanning_deletion_with_another_alternate_allele_is_rejected(files):
    with pytest.raises(PhaseSpaceError, match="spanning deletion combined"):
        load(files, [row(31, "CCC", "C"), row(32, "C", "*,T", gt="1/2")])


def test_multi_base_block_anchor_dp_does_not_establish_minimum_depth(files):
    locus = load(files, [row(1, "A", "<NON_REF>", gt="0/0", info="END=50")])
    assert all(item.status == "unknown" for item in locus.evidence)
    assert any("missing MIN_DP for multi-base reference block" in warning for warning in locus.warnings)


def test_multi_base_block_can_disable_minimum_depth_requirement(files):
    locus = load(files, [row(1, "A", "<NON_REF>", gt="0/0", info="END=50")], min_dp=0)
    assert all(item.status == "reference" for item in locus.evidence)


def test_structural_homref_end_is_not_a_reference_block(files):
    locus = load(files, [row(16, "A", "<DEL>", gt="0/0", info="END=50")])
    reference = [(item.start, item.end) for item in locus.evidence if item.status == "reference"]
    assert reference == [(15, 16)]


@pytest.mark.parametrize("thresholds", [{"min_gq": float("nan")}, {"min_gq": float("inf")}, {"min_dp": float("inf")}, {"min_dp": 1.5}, {"min_dp": -1}])
def test_quality_thresholds_must_be_finite_nonnegative(files, thresholds):
    with pytest.raises(PhaseSpaceError, match="finite and nonnegative"):
        load(files, [], **thresholds)


def test_declared_vcf_contig_mismatch_rejected(files):
    reference, annotation, write = files
    path = write([], header=HEADER.replace("ID=chr1", "ID=chr2"))
    with pytest.raises(PhaseSpaceError, match="absent from the VCF header"):
        load_locus(path, reference, annotation, "TX")


def test_declared_vcf_contig_length_mismatch_rejected(files):
    reference, annotation, write = files
    path = write([], header=HEADER.replace("length=90", "length=91"))
    with pytest.raises(PhaseSpaceError, match="Contig length mismatch"):
        load_locus(path, reference, annotation, "TX")


@pytest.mark.parametrize("translation_table", ["2", "11"])
def test_declared_nonstandard_translation_table_rejected(files, translation_table):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text().replace('gene_id "GENE";', f'transl_table "{translation_table}"; gene_id "GENE";'))
    with pytest.raises(PhaseSpaceError, match="unsupported translation table"):
        load_locus(write([]), reference, annotation, "TX")


def test_gff3_transcript_translation_table_rejected(files):
    reference, annotation, write = files
    annotation.write_text(
        "chr1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=TX;transl_table=2\n"
        "chr1\ttest\tCDS\t13\t18\t.\t+\t0\tParent=TX\n"
        "chr1\ttest\tCDS\t31\t39\t.\t+\t0\tParent=TX\n"
    )
    with pytest.raises(PhaseSpaceError, match="unsupported translation table"):
        load_locus(write([]), reference, annotation, "TX")


def test_gtf_transcript_translation_exception_rejected(files):
    reference, annotation, write = files
    annotation.write_text('chr1\ttest\ttranscript\t1\t90\t.\t+\t.\ttranscript_id "TX"; transl_except "(pos:31..33,aa:Sec)";\n' + annotation.read_text())
    with pytest.raises(PhaseSpaceError, match="transl_except"):
        load_locus(write([]), reference, annotation, "TX")


def test_gff3_cds_translation_exception_rejected(files):
    reference, annotation, write = files
    annotation.write_text(
        "chr1\ttest\tCDS\t13\t18\t.\t+\t0\tParent=TX\n"
        "chr1\ttest\tCDS\t31\t39\t.\t+\t0\tParent=TX;transl_except=%28pos:31..33,aa:Sec%29\n"
    )
    with pytest.raises(PhaseSpaceError, match="transl_except"):
        load_locus(write([]), reference, annotation, "TX")


def test_selected_selenocysteine_feature_rejected(files):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text() + 'chr1\ttest\tSelenocysteine\t31\t33\t.\t+\t.\ttranscript_id "TX";\n')
    with pytest.raises(PhaseSpaceError, match="Selenocysteine feature"):
        load_locus(write([]), reference, annotation, "TX")


def test_unselected_selenocysteine_feature_does_not_block_transcript(files):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text() + 'chr1\ttest\tSelenocysteine\t31\t33\t.\t+\t.\ttranscript_id "OTHER";\n')
    assert load_locus(write([]), reference, annotation, "TX").transcript.id == "TX"


def _custom(files, sequence, annotation_text):
    reference, annotation, write = files
    reference.write_text(">chr1\n" + sequence.ljust(90, "A") + "\n")
    Path(str(reference) + ".fai").unlink()
    pysam.faidx(str(reference))
    annotation.write_text(annotation_text)
    return reference, annotation, write([])


def gtf(feature, start, end, strand="+", transcript="TX", phase="0"):
    return f'chr1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t{phase}\tgene_id "G"; transcript_id "{transcript}";\n'


def test_gtf_stop_codon_is_joined_to_the_cds(files):
    sequence = "A" * 12 + "ATGAAA" + "A" * 12 + "CCCGGGTAA"
    reference, annotation, vcf = _custom(files, sequence, gtf("CDS", 13, 18) + gtf("CDS", 31, 36) + gtf("stop_codon", 37, 39))
    transcript = load_locus(vcf, reference, annotation, "TX").transcript
    assert transcript.cds_intervals == ((12, 18), (30, 39))
    assert transcript.reference_cds == "ATGAAACCCGGGTAA"


def test_gtf_stop_codon_split_across_exons(files):
    sequence = "A" * 10 + "ATGAAAT" + "C" * 10 + "AAG"
    reference, annotation, vcf = _custom(files, sequence, gtf("CDS", 11, 16) + gtf("stop_codon", 17, 17) + gtf("stop_codon", 28, 29, phase="2"))
    transcript = load_locus(vcf, reference, annotation, "TX").transcript
    assert transcript.cds_intervals == ((10, 17), (27, 29))
    assert transcript.reference_cds == "ATGAAATAA"


def test_reverse_strand_gtf_stop_codon_precedes_the_cds_genomically(files):
    sequence = "A" * 10 + "TTACCCCAT"
    reference, annotation, vcf = _custom(files, sequence, gtf("CDS", 14, 19, "-") + gtf("stop_codon", 11, 13, "-"))
    transcript = load_locus(vcf, reference, annotation, "TX").transcript
    assert transcript.cds_intervals == ((10, 19),)
    assert transcript.reference_cds == "ATGGGGTAA"


def test_stop_codon_already_inside_cds_is_not_added_twice(files):
    reference, annotation, write = files
    annotation.write_text(annotation.read_text() + gtf("stop_codon", 37, 39))
    assert load_locus(write([]), reference, annotation, "TX").transcript.reference_cds == "ATGAAACCCGGGTAA"


@pytest.mark.parametrize("requested", ["NM_000001.2", "NM_000001", "rna-NM_000001.2"])
def test_refseq_gff3_transcript_resolves_from_its_accession(files, requested):
    reference, annotation, write = files
    annotation.write_text(
        "##gff-version 3\n"
        "chr1\tRefSeq\tgene\t1\t90\t.\t+\t.\tID=gene-G;Name=G\n"
        "chr1\tRefSeq\tmRNA\t1\t90\t.\t+\t.\tID=rna-NM_000001.2;Parent=gene-G;transcript_id=NM_000001.2\n"
        "chr1\tRefSeq\tCDS\t13\t18\t.\t+\t0\tID=cds-NP_1;Parent=rna-NM_000001.2\n"
        "chr1\tRefSeq\tCDS\t31\t39\t.\t+\t0\tID=cds-NP_1;Parent=rna-NM_000001.2\n"
    )
    locus = load_locus(write([]), reference, annotation, requested)
    assert locus.transcript.gene == "G"
    assert locus.transcript.reference_cds == "ATGAAACCCGGGTAA"


def test_ensembl_gff3_versioned_id_uses_the_version_attribute(files):
    reference, annotation, write = files
    annotation.write_text(
        "chr1\tens\tmRNA\t1\t90\t.\t+\t.\tID=transcript:ENST1;transcript_id=ENST1;version=4\n"
        "chr1\tens\tCDS\t13\t18\t.\t+\t0\tParent=transcript:ENST1\n"
        "chr1\tens\tCDS\t31\t39\t.\t+\t0\tParent=transcript:ENST1\n"
    )
    for requested in ("ENST1", "ENST1.4", "transcript:ENST1"):
        assert load_locus(write([]), reference, annotation, requested).transcript.reference_cds == "ATGAAACCCGGGTAA"
    with pytest.raises(PhaseSpaceError, match="No CDS features"):
        load_locus(write([]), reference, annotation, "ENST1.3")


def test_unversioned_id_matching_several_versions_is_rejected(files):
    reference, annotation, write = files
    annotation.write_text(
        'chr1\ttest\ttranscript\t1\t90\t.\t+\t.\ttranscript_id "TX.1";\n'
        'chr1\ttest\ttranscript\t1\t90\t.\t+\t.\ttranscript_id "TX.2";\n'
        + gtf("CDS", 13, 18, transcript="TX.1") + gtf("CDS", 31, 39, transcript="TX.1")
        + gtf("CDS", 13, 18, transcript="TX.2") + gtf("CDS", 31, 39, transcript="TX.2")
    )
    with pytest.raises(PhaseSpaceError, match="several annotated versions"):
        load_locus(write([]), reference, annotation, "TX")
    assert load_locus(write([]), reference, annotation, "TX.2").transcript.id == "TX.2"


def phase_file(tmp_path, records, samples=("SAMPLE",)):
    path = tmp_path / "phased.vcf"
    header = HEADER.replace("\tSAMPLE\n", "\t" + "\t".join(samples) + "\n")
    path.write_text(header + "".join(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:PS\t{values}\n" for pos, ref, alt, values in records))
    return path


def test_phase_from_supplies_gt_order_and_phase_set(files, tmp_path):
    reference, annotation, write = files
    phased = phase_file(tmp_path, [(16, "A", "T", "1|0:5"), (32, "C", "T", "0|1:5")])
    locus = load_locus(write([row(16, "A", "T"), row(32, "C", "T")]), reference, annotation, "TX", phase_from=phased)
    first, second = locus.variants
    assert (first.gt, first.alleles, first.cds_alleles, first.phased, first.phase_set) == ((1, 0), ("T", "A"), ("T", "A"), True, "5")
    assert first.id == "chr1:16:A>T"
    assert (second.gt, second.alleles, second.phased, second.phase_set) == ((0, 1), ("C", "T"), True, "5")
    assert any("2 of 2 heterozygous coding calls phased" in w for w in locus.warnings)
    assert locus.source == "calls.vcf + phase from phased.vcf"
    # Opposite orientation within one phase set leaves a single arrangement.
    from phasespace.analysis import analyze
    result = analyze(locus)
    assert result["scenario_count"] == 1
    assert sorted(len(h["variant_ids"]) for h in result["scenarios"][0]["haplotypes"]) == [1, 1]


def test_phase_from_replaces_the_inputs_own_phase(files, tmp_path):
    reference, annotation, write = files
    phased = phase_file(tmp_path, [(16, "A", "T", "0|1:5")])
    calls = [row(16, "A", "T", fmt="GT:GQ:DP:PS", values="1|0:50:12:9"), row(32, "C", "T", fmt="GT:GQ:DP:PS", values="0|1:50:12:9")]
    first, second = load_locus(write(calls), reference, annotation, "TX", phase_from=phased).variants
    assert (first.gt, first.phase_set) == ((0, 1), "5")
    assert (second.phased, second.phase_set) == (False, None)


def test_phase_from_matches_alleles_by_sequence_not_index(files, tmp_path):
    reference, annotation, write = files
    phased = phase_file(tmp_path, [(16, "A", "G,T", "1|0:5")])
    variant = load_locus(write([row(16, "A", "T,G", gt="0/2")]), reference, annotation, "TX", phase_from=phased).variants[0]
    assert (variant.gt, variant.alleles, variant.phased) == ((2, 0), ("G", "A"), True)


def test_phase_from_ignores_records_with_different_alleles(files, tmp_path):
    reference, annotation, write = files
    phased = phase_file(tmp_path, [(16, "A", "C", "0|1:5")])
    locus = load_locus(write([row(16, "A", "T")]), reference, annotation, "TX", phase_from=phased)
    assert not locus.variants[0].phased
    assert any("0 of 1 heterozygous coding call phased; 1 not phased there is treated as unphased" in w for w in locus.warnings)


def test_phase_from_sample_selection(files, tmp_path):
    reference, annotation, write = files
    single = phase_file(tmp_path, [(16, "A", "T", "0|1:5")], samples=("RENAMED",))
    locus = load_locus(write([row(16, "A", "T")]), reference, annotation, "TX", phase_from=single)
    assert locus.variants[0].phased and any("sample RENAMED" in w for w in locus.warnings)
    several = phase_file(tmp_path, [(16, "A", "T", "0|1:5\t0|1:5")], samples=("A", "B"))
    with pytest.raises(PhaseSpaceError, match="absent from the phase file"):
        load_locus(write([row(16, "A", "T")]), reference, annotation, "TX", phase_from=several)


def test_vcf_contig_without_chr_prefix_is_matched(files):
    reference, annotation, write = files
    path = write([row(16, "A", "T").replace("chr1\t", "1\t", 1)], header=HEADER.replace("ID=chr1", "ID=1"))
    locus = load_locus(path, reference, annotation, "TX")
    assert locus.transcript.contig == "1"
    assert locus.variants[0].id == "1:16:A>T"
    assert any("Contig names differ" in w and "VCF '1'" in w for w in locus.warnings)


def test_reference_contig_without_chr_prefix_is_matched(files):
    reference, annotation, write = files
    reference.write_text(reference.read_text().replace(">chr1", ">1"))
    Path(str(reference) + ".fai").unlink()
    pysam.faidx(str(reference))
    # The VCF header names chr1; with the FASTA as "1", lengths are still compared.
    locus = load_locus(write([row(16, "A", "T")]), reference, annotation, "TX")
    assert locus.transcript.reference_cds == "ATGAAACCCGGGTAA"
    assert locus.variants[0].id == "chr1:16:A>T"


def test_contig_aliases_cover_chr_prefix_and_mitochondria():
    from phasespace.inputs import _contig_aliases
    assert _contig_aliases("chr7") == ["chr7", "7"]
    assert _contig_aliases("X") == ["X", "chrX"]
    assert set(_contig_aliases("MT")) == {"MT", "chrMT", "chrM", "M"}


def test_variant_id_lists_alternate_alleles_only(files):
    reference, annotation, write = files
    multi = load(files, [row(16, "A", "T,G", gt="2|1", fmt="GT:GQ:DP:PS", values="2|1:50:12:3")]).variants[0]
    assert multi.id == "chr1:16:A>G,T"
    assert multi.alleles == ("G", "T")


def test_refcall_reference_genotype_is_reference_evidence(files):
    calls = [row(16, "A", "T", gt="0/0", filt="RefCall"), row(32, "C", "T", gt="0/0", filt="RefCall", values="0/0:5:12")]
    locus = load(files, calls)
    assert [(e.start, e.end, e.status) for e in locus.evidence if e.start in (15, 31)] == [(15, 16, "reference"), (31, 32, "unknown")]
    assert any("GQ below 20" in w for w in locus.warnings)
    assert not locus.variants


def test_refcall_on_a_variant_genotype_still_excludes_it(files):
    locus = load(files, [row(16, "A", "T", gt="0/1", filt="RefCall")])
    assert not locus.variants
    assert any("FILTER=RefCall" in w for w in locus.warnings)


def test_splice_region_calls_are_flagged_not_modelled(files):
    # Coding exons [12,18) and [30,39): junctions at 18 and 30.
    calls = [row(5, "A", "G"), row(17, "A", "G"), row(21, "A", "G"), row(28, "A", "C", values="0/1:5:12"), row(29, "A", "T", gt="1/1")]
    locus = load(files, calls)
    hits = {hit["id"]: (hit["where"], hit["distance"], hit["junction"]) for hit in locus.splice_region}
    assert hits == {"chr1:17:A>G": ("exon", 2, 18), "chr1:21:A>G": ("intron", 3, 18), "chr1:29:A>T": ("intron", 2, 30)}
    assert [v.id for v in locus.variants] == ["chr1:17:A>G"]  # Coding consequence still modelled.
    assert any("3 called variants lie in a splice region" in w for w in locus.warnings)
    from phasespace.analysis import analyze
    result = analyze(locus)
    assert len(result["splice_region"]) == 3 and result["summary"]["splice_region"] == 3


def test_single_exon_transcript_has_no_splice_region(files):
    reference, annotation, write = files
    annotation.write_text('chr1\ttest\tCDS\t31\t39\t.\t+\t0\tgene_id "G"; transcript_id "TX";\n')
    reference.write_text(reference.read_text().replace("CCCGGGTAA", "ATGGGGTAA"))
    Path(str(reference) + ".fai").unlink()
    pysam.faidx(str(reference))
    assert load_locus(write([row(29, "A", "T"), row(31, "A", "G")]), reference, annotation, "TX").splice_region == []
