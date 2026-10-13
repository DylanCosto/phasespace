# PhaseSpace

**Explore how uncertain phasing changes the possible coding consequences of your variants.**

PhaseSpace reads a local VCF/gVCF or BCF, a matching indexed reference FASTA, and a GTF/GFF3 transcript annotation. It generates a standalone interactive HTML report for one sample and one transcript. No server, account, cloud upload, or external JavaScript is required to view the report.

Two heterozygous stop variants can occur on the same gene copy or different copies. PhaseSpace retains both arrangements when the data do not resolve that relationship. The viewer lets you explore hypothetical phase constraints, examine the two reconstructed coding sequences, and identify which block relationships distinguish different outcomes.

## Quick start

Python 3.10+ on Linux or macOS:

```bash
git clone https://github.com/DylanCosto/phasespace.git
cd phasespace
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
phasespace demo --out reports/demo.html --json reports/demo.json
```

Open `reports/demo.html` in a browser. The demo uses completely synthetic sequences and includes an intentional gap in reference-call evidence. To keep the underlying example inputs:

```bash
phasespace demo --out reports/demo.html --force --write-inputs demo-inputs
```

## Analyze your data

```bash
phasespace report \
  --vcf sample.g.vcf.gz \
  --fasta reference.fa \
  --annotation genes.gtf \
  --transcript ENST00000000000.1 \
  --sample SAMPLE_ID \
  --out reports/transcript.html \
  --json reports/transcript.json
```

Use the **exact transcript ID** in the annotation. The reference build and contig names must match across all inputs. A `.fai` FASTA index is required; prepare it with `samtools faidx reference.fa`. A bgzip VCF with a tabix/CSI index, or indexed BCF, enables regional reading. Unindexed VCFs are streamed with an explicit notice. The annotation is scanned to locate the requested transcript.

Default call thresholds are GQ ≥ 20 and depth ≥ 5. Multi-base gVCF reference blocks require `MIN_DP` for a positive depth threshold: anchor `DP` alone cannot establish depth throughout a block. Single-base records use `DP` (or `MIN_DP` for a reference block). Missing required quality values, missing genotypes, and failed filters produce unknown evidence. Set thresholds with `--min-gq` and `--min-dp`. Setting minimum depth to zero permits a block with only `DP`, explicitly waiving the block-wide depth check; missing quality/depth fields still remain unknown. Multi-sample files require `--sample`.

## What the report shows

- All supported phase arrangements consistent with the recorded genotype order and phase sets.
- Distinct **unordered pairs of coding sequences**, with the original arrangements retained.
- Joint sequence effects, including premature stops, amino-acid changes, frameshifts, and compensating indels.
- Reference-call support, called-variant evidence, and unknown intervals across the annotated CDS.
- Hypothetical block-orientation constraints and informative block-pair comparisons.
- A local JSON export with the analyzed sequences, calls, evidence, and scenarios.

Hypothetical constraints do not change the VCF or become observed evidence. Each scenario is a possible consequence of the **included called variants**. Unknown intervals permit additional unobserved variants, so these are not exhaustive possible biological sequences. “Reference-identical CDS” counts exact sequence matches within this model; it does not establish protein function.

## Scope of version 0.1

This is a research prototype for diploid coding loci, using the standard genetic code. It supports SNVs, substitutions and sequence-resolved small indels contained within a single CDS segment, including multi-allelic diploid genotypes. Calls on the negative strand are projected into transcript orientation before applying edits.

The first CDS phase must be zero, subsequent CDS phases must agree, and joined CDS length must be divisible by three. Only `CDS` annotation features are used. In GTF files that exclude the terminal stop codon from CDS, terminal stop-loss assessment is consequently unavailable. Partial models, symbolic called alleles, spanning deletions, overlapping edits, and edits crossing CDS boundaries are rejected rather than approximated. Reference mismatches are errors.

The default cap is eight independent heterozygous blocks (at most 128 arrangements after exchanging the two copy labels). `--max-blocks` can raise this to twelve. Larger problems stop with an error; scenarios are never silently sampled or truncated. The report has no scenario probabilities. Phased genotypes without `PS` share an implicit phase set as specified by VCF; different explicit phase sets remain independently oriented.

Splicing, structural variants, copy-number changes, mitochondrial genetic codes, recoded amino acids such as selenocysteine, protein-domain annotation, nonsense-mediated decay, expression, and clinical pathogenicity are outside this version's model. Explicit nonstandard translation annotations are rejected. The tool does not perform phasing or establish a diagnosis.

## Privacy and reproducibility

Analysis is local. The HTML contains the selected sample identifier, locus calls, reconstructed sequences, and evidence, so treat it like the source genomic data when sharing. There is no telemetry and the report makes no network requests. Files in `examples/` are synthetic. The JSON includes software version and analysis duration.

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest
python -m build
```

See [the model and validation notes](docs/model.md) for details and [the validation record](docs/validation.md) for the checks performed on this first version.

## Related work

[GEMPROT](https://pubmed.ncbi.nlm.nih.gov/30508040/) and [ProHap Explorer](https://github.com/ProGenNo/ProHapExplorer) visualize haplotype-level protein variation. [BCFtools/csq](https://samtools.github.io/bcftools/howtos/csq-calling.html) computes haplotype-aware consequences. PhaseSpace's focus is interactive exploration of unresolved arrangements together with explicit reference-call evidence. This is a workflow differentiation, not a claim that haplotype-aware consequence analysis is new.

Input semantics follow the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf); indexed file access uses [pysam](https://pysam.readthedocs.io/).
