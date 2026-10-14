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

Give the transcript ID as written in the annotation, with or without its version (for example `ENST00000371998.8` or `ENST00000371998`). GFF3 feature IDs such as `transcript:ENST…` and RefSeq `rna-NM_…` also resolve, as does a RefSeq accession on its own. An unversioned ID that matches more than one annotated version is rejected. All inputs must use the same reference build; contig names may differ only by a `chr` prefix. A `.fai` FASTA index is required; prepare it with `samtools faidx reference.fa`. A bgzip VCF with a tabix/CSI index, or indexed BCF, enables regional reading. Unindexed VCFs are streamed with an explicit notice. The annotation is scanned to locate the requested transcript.

Default call thresholds are GQ ≥ 20 and depth ≥ 5. Multi-base gVCF reference blocks require `MIN_DP` for a positive depth threshold: anchor `DP` alone cannot establish depth throughout a block. Single-base records use `DP` (or `MIN_DP` for a reference block). Missing required quality values, missing genotypes, and failed filters produce unknown evidence. Set thresholds with `--min-gq` and `--min-dp`. Setting minimum depth to zero permits a block with only `DP`, explicitly waiving the block-wide depth check; missing quality/depth fields still remain unknown. Multi-sample files require `--sample`. DeepVariant `RefCall` records with a `0/0` genotype count as reference evidence when they pass these thresholds. Contig names may differ between inputs by a `chr` prefix; positions in the report use the VCF's names.

Deletions in gVCFs are handled whichever way the writer records the deleted bases. A reference block that continues through the bases removed by a called deletion is not treated as a conflict, and a `*` spanning-deletion allele is accepted when an upstream called deletion covers it. A `*` allele combined with another alternate allele at the same site is rejected, because placing it requires read-level phasing.

### Pick a transcript by gene

Use `--gene SYMBOL` (or a gene ID) instead of `--transcript`. PhaseSpace picks the gene's MANE Select transcript, then its Ensembl canonical transcript, then its RefSeq Select transcript, and says which it chose. If none is tagged, it lists the transcripts so you can pick one with `--transcript`.

### Show annotations you already have

`--info CLNSIG,AF` shows those INFO fields for each variant in reports and scans. By default they come from the main VCF. `--info-from clinvar.vcf.gz` reads them from a separate indexed VCF, such as ClinVar or gnomAD sites, matching calls by position, REF and ALT. Contig names may differ by a `chr` prefix. Per-allele fields (`Number=A` or `R`) are taken for the called alternate alleles. Nothing is looked up online.

### Scan every gene

```bash
phasespace scan --vcf sample.g.vcf.gz --fasta reference.fa --annotation genes.gtf.gz \
  --sample SAMPLE_ID --out-dir scan/
```

This analyses one representative transcript per protein-coding gene and ranks the genes:

1. **Phase decides** whether a copy keeps the reference protein.
2. **No reference protein, disruptive:** no arrangement keeps a copy with the reference protein, and a frameshift, premature stop, or start or stop loss is involved.
3. **No reference protein** (often common homozygous missense variants).
4. **Reference protein kept.**

The scan writes `scan/index.html` (sortable and filterable), `scan.tsv`, `scan.json`, and reports for tiers 1–2 (`--reports all` writes one for every gene with variants). Restrict it to a panel with `--genes CFTR,PAH` or `--genes @panel.txt`. `--phase-from`, thresholds and `--max-blocks` work as for `report`. Genes on contigs absent from the reference or VCF are counted but not listed, so single-chromosome and exome inputs work. A whole-genome GENCODE scan of HG002 chr20 calls takes under 20 seconds on 12 processes.

### Using phase from WhatsHap, longphase or a trio

Phasing tools write variant-only VCFs. In such a file, absent records are unknown evidence, never reference. Pass the phased VCF separately so the gVCF still supplies the evidence:

```bash
phasespace report --vcf sample.g.vcf.gz --phase-from phased.vcf.gz \
  --fasta reference.fa --annotation genes.gtf --transcript ENST00000000000.1 \
  --sample SAMPLE_ID --out reports/transcript.html
```

Heterozygous calls are matched by position, REF and called allele sequences, so allele numbering may differ between the files (for example in a multi-sample trio VCF). A matched call takes its GT order and `PS` from the phase file. With `--phase-from`, phase comes only from that file: any heterozygous call it does not phase is treated as unphased, so phase sets from the two files are never mixed. The report notes how many calls were phased. If the phase file has one sample with a different name, that sample is used. `scripts/transfer_phase.py` writes the same phase into a gVCF, for use with other tools.

## What the report shows

- All supported phase arrangements consistent with the recorded genotype order and phase sets.
- Distinct **unordered pairs of coding sequences**, with the original arrangements retained.
- Joint sequence effects, including premature stops, amino-acid changes, frameshifts, and compensating indels.
- Reference-call support, called-variant evidence, and unknown intervals across the annotated CDS.
- Hypothetical block-orientation constraints and informative block-pair comparisons.
- A verdict on whether any arrangement keeps a copy with the reference protein sequence, and whether that depends on phase. Synonymous changes alter the CDS but not this verdict.
- Both gene copies side by side, with each copy's alternate alleles, its predicted protein length on a shared scale, and the first stop.
- Per-block flip controls that move a block's alleles to the other copy.
- A local JSON export with the analyzed sequences, calls, evidence, and scenarios.

The report stores each copy as its called edits rather than as full sequences, and rebuilds sequences in the browser on demand. A 10 kb CDS at the 12-block cap produces a report of about 1.6 MB. The report follows the system light or dark setting.

Hypothetical constraints do not change the VCF or become observed evidence. Each scenario is a possible consequence of the **included called variants**. Unknown intervals permit additional unobserved variants, so these are not exhaustive possible biological sequences. “Reference-identical CDS” counts exact sequence matches within this model; it does not establish protein function.

## Scope of version 0.1

This is a research prototype for diploid coding loci, using the standard genetic code. It supports SNVs, substitutions and sequence-resolved small indels contained within a single CDS segment, including multi-allelic diploid genotypes. Calls on the negative strand are projected into transcript orientation before applying edits.

The first CDS phase must be zero, subsequent CDS phases must agree, and joined CDS length must be divisible by three. `CDS` features define the coding sequence. GTF files exclude the terminal stop codon from `CDS`, so their `stop_codon` features are joined to it, including stop codons split across exons. Partial models, symbolic called alleles, spanning-deletion alleles not covered by a called upstream deletion, overlapping edits, and edits crossing CDS boundaries are rejected rather than approximated. Reference mismatches are errors.

The default cap is eight independent heterozygous blocks (at most 128 arrangements after exchanging the two copy labels). `--max-blocks` can raise this to twelve. Larger problems stop with an error; scenarios are never silently sampled or truncated. The report has no scenario probabilities. Phased genotypes without `PS` share an implicit phase set as specified by VCF; different explicit phase sets remain independently oriented.

Splicing effects are not modelled. Instead, passing calls within 3 exonic or 8 intronic bases of a junction between coding exons are flagged as splice-region variants in reports and scans, so a reference-protein verdict is never read as clearing them. Structural variants, copy-number changes, mitochondrial genetic codes, recoded amino acids such as selenocysteine, protein-domain annotation, nonsense-mediated decay, expression, and clinical pathogenicity are outside this version's model. Explicit nonstandard translation annotations are rejected. The tool does not perform phasing or establish a diagnosis.

## Privacy and reproducibility

Analysis is local. The HTML contains the selected sample identifier, locus calls, reconstructed sequences, and evidence, so treat it like the source genomic data when sharing. There is no telemetry and the report makes no network requests. Files in `examples/` are synthetic. The JSON includes software version and analysis duration.

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest
python -m build
```

See [the architecture guide](docs/architecture.md) for how the code is organised, [the model notes](docs/model.md) for interpretation, and [the validation record](docs/validation.md) for the checks performed, including on real HG002 data.

## Related work

[GEMPROT](https://pubmed.ncbi.nlm.nih.gov/30508040/) and [ProHap Explorer](https://github.com/ProGenNo/ProHapExplorer) visualize haplotype-level protein variation. [BCFtools/csq](https://samtools.github.io/bcftools/howtos/csq-calling.html) computes haplotype-aware consequences. PhaseSpace's focus is interactive exploration of unresolved arrangements together with explicit reference-call evidence. This is a workflow differentiation, not a claim that haplotype-aware consequence analysis is new.

Input semantics follow the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf); indexed file access uses [pysam](https://pysam.readthedocs.io/).
