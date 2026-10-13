# Synthetic example

All sequences, identifiers and calls in this directory were generated for PhaseSpace. They are not patient, participant, or organism data.

`synthetic/reference.fa`, `transcripts.gtf`, `sample.vcf`, and `sample.g.vcf` describe a 96-base CDS in three exons. Two unphased heterozygous substitutions introduce predicted stop codons at amino-acid positions 8 and 21. The gVCF has a six-base interval without evidence; the variant-only VCF has no reference blocks.

Generate an indexed copy and the interactive report from an installed package:

```bash
phasespace demo --out reports/demo.html --write-inputs demo-inputs
```

Or analyze these tracked plain-text inputs directly:

```bash
phasespace report --vcf examples/synthetic/sample.g.vcf \
  --fasta examples/synthetic/reference.fa \
  --annotation examples/synthetic/transcripts.gtf \
  --transcript DEMO1.1 --out reports/example.html
```

The second command streams the small unindexed gVCF and reports that fact. It produces the same two phase arrangements as the indexed demo.
