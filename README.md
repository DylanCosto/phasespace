# PhaseSpace

[![tests](https://github.com/DylanCosto/phasespace/actions/workflows/tests.yml/badge.svg)](https://github.com/DylanCosto/phasespace/actions/workflows/tests.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Explore how uncertain phasing changes the possible coding consequences of your variants.**

Two heterozygous loss-of-function variants on the same gene copy leave the other copy intact; on different copies they leave none. When the data don't resolve phase, PhaseSpace keeps every arrangement consistent with the calls, rebuilds and translates both copies jointly, and tells you whether a copy keeps the reference protein in every arrangement, in none, or only depending on phase. If reads can settle it, it shows which reads.

It works from local files only: a VCF/gVCF or BCF, an indexed reference FASTA and a GTF/GFF3 annotation, plus optional phased VCF and BAM/CRAM files. Results are self-contained interactive HTML reports that make no network requests, with JSON and TSV for pipelines. Run it on one gene (`report`) or every gene (`scan`).

## Quick start

Python 3.10+ on Linux or macOS:

```bash
python -m pip install phasespace-vcf          # after the 0.2.0 release; or install from a clone:
git clone https://github.com/DylanCosto/phasespace.git && python -m pip install ./phasespace
phasespace demo --out reports/demo.html --json reports/demo.json
```

A Docker image (`ghcr.io/dylancosto/phasespace`) and Nextflow modules are also available; see [Docker and Nextflow](#docker-and-nextflow).

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

DeepVariant 1.10 and GATK 4 HaplotypeCaller gVCFs (both checked on HG002) start the next reference block after a deletion's REF span. For writers that instead continue a reference block through deleted bases, that block is not treated as a conflict. A `*` spanning-deletion allele, as in jointly genotyped files, is accepted when a called upstream deletion covers it. A `*` allele combined with another alternate allele at the same site is rejected, because placing it requires read-level phasing.

### Pick a transcript by gene

Use `--gene SYMBOL` (or a gene ID) instead of `--transcript`. PhaseSpace picks the gene's MANE Select transcript, then its Ensembl canonical transcript, then its RefSeq Select transcript, and says which it chose. If none is tagged, it lists the transcripts so you can pick one with `--transcript`.

### Read evidence from BAM/CRAM

Add `--reads sample.bam` (repeatable; CRAM uses `--fasta`) to `report` or `scan`. For every read or read pair that covers two heterozygous calls in different blocks, PhaseSpace records which alleles it carries. That shows whether the blocks share an orientation (alternate alleles on one copy) or not. Each fragment counts once per block pair. SNVs use the aligned base (base quality ≥ 13), and indels compare the read with both expected sequences across a flanked window. Reads need MAPQ ≥ 20 and must not be secondary, duplicate or QC-fail.

The report shows these as **observed counts**, kept separate from hypotheses. A relationship is "supported" with at least 3 fragments and 90% agreement, and **Apply read-supported relationships** adds them as hypotheses tagged "reads". Within existing phase blocks, the report also shows how often reads agree with the input phase. Scans add `verdict_with_reads`. Nothing is converted into a probability, and the input calls never change. Long reads matter most: on HG002, HiFi reads resolved CDC25B variants 1–6 kb apart that Illumina pairs could not reach.

### Inheritance evidence from parents

`--parents FATHER,MOTHER` (sample names; `--parents-from trio.vcf.gz` if they are not in the main VCF) works out which parent passed on each allele of the child's heterozygous calls. Calls whose alleles came from the same parent sit on one copy. A site is informative when only one assignment fits the parents' genotypes. Sites where both parents are heterozygous stay open, and impossible combinations are counted as Mendelian conflicts and never used. Parents' genotypes need passing filters and GQ ≥ `--min-gq`. Depth is checked when present; merged trio files often lack it for genotypes filled from gVCF blocks. Missing parental sites can come from gVCF reference blocks. As with reads, the report shows this as observed evidence with an **Apply inherited relationships** button, and scans add `verdict_with_parents` and `mendelian_conflicts`.

### Link blocks by evidence before enumerating (`--link-evidence`)

Without it, read and inheritance evidence only filter arrangements after they are enumerated, so a gene over the block cap cannot be analysed at all. `--link-evidence` (with `--reads` and/or `--parents`) merges blocks the evidence clearly links before enumerating: the same "supported" thresholds, strongest evidence first. A relationship that contradicts stronger ones is skipped and reported. Linked blocks are labelled with their sources, and the report states that the arrangements assume this evidence. On HG002, SIRPA (20 unphased heterozygous calls, far over the cap) becomes a single arrangement, and HiFi reads, Illumina reads and the parents each give the identical result.

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

### Docker and Nextflow

```bash
docker run --rm -v "$PWD:/data" -w /data ghcr.io/dylancosto/phasespace:0.2.0 \
  phasespace scan --vcf sample.g.vcf.gz --fasta ref.fa --annotation genes.gtf.gz --out-dir scan
```

The image (see `Dockerfile`) is published to GitHub Container Registry for each release. `nextflow/` provides nf-core-style `PHASESPACE_SCAN` and `PHASESPACE_REPORT` modules and an example samplesheet pipeline (`nextflow run nextflow/main.nf -profile docker,test`); see [nextflow/README.md](nextflow/README.md).

## What the report shows

- All supported phase arrangements consistent with the recorded genotype order and phase sets.
- Distinct **unordered pairs of coding sequences**, with the original arrangements retained.
- Joint sequence effects, including premature stops, amino-acid changes, frameshifts, and compensating indels.
- HGVS-style protein changes for each copy (`p.Gln8Ter`, `p.Val7_Leu8del`, `p.Lys2AsnfsTer2`), derived from the copy's jointly edited sequence. Two calls in one codon give the amino acid actually made. On HG002, SIRPA's c.298A>G (alone p.Asn100Asp) and c.300C>A (alone p.Asn100Lys) sit on one copy and give p.Asn100Glu. Substitutions follow HGVS exactly; indels are trimmed but not shifted 3′ across repeats outside the changed stretch.
- Reference-call support, called-variant evidence, and unknown intervals across the annotated CDS.
- Hypothetical block-orientation constraints and informative block-pair comparisons.
- A verdict on whether any arrangement keeps a copy with the reference protein sequence, and whether that depends on phase. Synonymous changes alter the CDS but not this verdict.
- Both gene copies side by side, with each copy's alternate alleles, its predicted protein length on a shared scale, and the first stop.
- Per-block flip controls that move a block's alleles to the other copy.
- A local JSON export with the analyzed sequences, calls, evidence, and scenarios.

The report stores each copy as its called edits rather than as full sequences, and rebuilds sequences in the browser on demand. A 10 kb CDS at the 12-block cap produces a report of about 1.6 MB. The report follows the system light or dark setting.

Hypothetical constraints do not change the VCF or become observed evidence. Each scenario is a possible consequence of the **included called variants**. Unknown intervals permit additional unobserved variants, so these are not exhaustive possible biological sequences. “Reference-identical CDS” counts exact sequence matches within this model; it does not establish protein function.

## Scope

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

## Citing

If you use PhaseSpace, please cite it using [CITATION.cff](CITATION.cff) (GitHub shows a "Cite this repository" button). A JOSS paper is in preparation (`paper/`).

## Related work

[GEMPROT](https://pubmed.ncbi.nlm.nih.gov/30508040/) and [ProHap Explorer](https://github.com/ProGenNo/ProHapExplorer) visualize haplotype-level protein variation. [BCFtools/csq](https://samtools.github.io/bcftools/howtos/csq-calling.html) computes haplotype-aware consequences. PhaseSpace's focus is interactive exploration of unresolved arrangements together with explicit reference-call evidence. This is a workflow differentiation, not a claim that haplotype-aware consequence analysis is new.

Input semantics follow the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf); indexed file access uses [pysam](https://pysam.readthedocs.io/).
