# HG002 chr20 regression fixture

Real GIAB HG002 calls for six GRCh38 chr20 transcripts, cut by `scripts/make_hg002_fixture.py`.
Coordinates are unchanged; reference bases outside the coding spans (plus 100 bp) are masked to N.

| Transcript | Why it is here |
|---|---|
| ENST00000371998.8 | NCOA3: heterozygous in-frame deletion, forward strand |
| ENST00000372726.5 | TP53TG5: homozygous in-frame deletion, reverse strand |
| ENST00000245960.10 | CDC25B: three unphased synonymous SNVs; CDS changes but protein never does |
| ENST00000217428.7 | SPINT3: three unphased SNVs including missense, reverse strand |
| ENST00000382398.4 | DEFB126: homozygous frameshift that reads through the stop codon |
| ENST00000358771.5 | SIRPA: 20 heterozygous variants; exceeds the block cap unless trio-phased |

Sources (all public):

- `pacbio.HG002.vcf.gz`: DeepVariant 1.10 DeepTrio PacBio chr20 gVCF, `gs://deepvariant/case-study-outputs/1.10.0/deeptrio/pacbio-chr20/HG002.child.g.vcf.gz`.
- `deepvariant.HG002.vcf.gz`: DeepVariant 1.10 DeepTrio gVCF published under `case-study-outputs/1.10.0/deeptrio/ont-r104-chr20/HG002.child.g.vcf.gz`. Despite the folder name, the file also holds exome-scale calls on other chromosomes, so its sequencing platform is not asserted here.
- `giab.HG002.vcf.gz`: GIAB HG002 NISTv4.2.1 GRCh38 benchmark VCF (unphased, variant-only).
- `trio.HG002.vcf.gz`: the `deepvariant` gVCF with phase from transmission. The DeepTrio HG002, HG003 and HG004 gVCFs from the same run were merged (`bcftools merge --gvcf`), genotypes with GQ < 20 set missing, phased with WhatsHap 2.8 `phase --ped` without reads, and copied back by `scripts/transfer_phase.py`. Sites heterozygous in all three samples stay unphased.
- `whatshap.HG002.vcf.gz`: the WhatsHap output itself (variant-only, all three samples), for `--phase-from`.
- `GRCh38.chr20.masked.fa.gz`: chr20 of the GRCh38 no-alt analysis set used by DeepVariant.
- `gencode.v46.gtf`: GENCODE v46 comprehensive annotation. `ensembl112.gff3`: Ensembl 112 chr20 GFF3, contig renamed from `20` to `chr20`.
