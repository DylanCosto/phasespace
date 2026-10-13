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
    path = write([], header=HEADER.replace("ID=chr1", "ID=1"))
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
