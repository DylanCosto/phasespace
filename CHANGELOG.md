# Changelog

## 0.2.0 — 2026-10-14

### Added
- `phasespace scan`: one representative transcript per protein-coding gene (MANE Select, then Ensembl canonical, then RefSeq Select), ranked by whether phase decides if a copy keeps the reference protein. Writes `scan.tsv`, `scan.json`, a sortable `index.html` and reports for candidates. Accepts `--genes` lists or files and runs in parallel.
- `--gene` for `report`, choosing the same representative transcript.
- `--reads BAM/CRAM`: read-backed evidence. Counts fragments supporting the same or opposite orientation for each block pair, and checks the input phase within blocks. Shown as observed counts, with an explicit apply button.
- `--phase-from`: take phase from a separate phased VCF (WhatsHap, longphase or trio) while the gVCF keeps the evidence.
- `--info` and `--info-from`: show INFO annotations from the VCF or a separate local VCF such as ClinVar or gnomAD sites.
- Splice-region flags for passing calls near coding-exon junctions, which are not modelled.
- Docker image, nf-core-style Nextflow modules with an example pipeline, MIT license.
- Real-data regression fixture (GIAB HG002 chr20) and validation on DeepVariant, GATK and GIAB inputs, HiFi and Illumina reads, and trio phasing.

### Changed
- The redesigned report leads with a protein-level verdict and shows the two copies side by side. Large analyses embed edits instead of sequences: 1.6 MB instead of 44 MB at the block cap.
- GTF `stop_codon` features are joined to the CDS, so stop loss is assessed for GENCODE and Ensembl GTFs.
- Transcript IDs resolve with or without version, and as GFF3/RefSeq IDs. Contig names may differ by a `chr` prefix.
- DeepVariant `RefCall` `0/0` calls count as reference evidence when they pass thresholds.
- Variant IDs list sorted alternate alleles only.
- Faster: engine ~3× faster, and annotation lookup ~14× faster (59 s to 4.2 s on GENCODE).

### Fixed
- gVCF deletions with overlapping reference blocks or `*` alleles no longer fail.
- A frameshift that gains a premature stop is no longer also labelled `stop_lost`.

## 0.1.0 — 2026-10-12

First version: one transcript, local VCF/gVCF, interactive standalone report.
