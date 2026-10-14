---
title: 'PhaseSpace: exploring how unresolved phase changes the coding consequences of called variants'
tags:
  - Python
  - genomics
  - phasing
  - compound heterozygosity
  - long reads
  - variant interpretation
authors:
  - name: Dylan Costo
    orcid: 0000-0000-0000-0000  # TODO: replace with your ORCID
    affiliation: 1
affiliations:
  - name: TODO affiliation
    index: 1
date: 14 October 2026
bibliography: paper.bib
---

# Summary

A person carries two copies of each autosomal gene. When a gene contains two or more heterozygous variants, their effect on the protein depends on phase: whether the variants sit on the same copy or on different copies. Two loss-of-function variants on one copy leave the other copy intact; on different copies they leave none. Variant callers often leave phase unresolved, and consequence predictors usually assess variants one at a time [@McLaren2016] or assume a phase [@Danecek2017].

PhaseSpace makes this uncertainty explicit. For a transcript, it enumerates every arrangement of the called heterozygous variants that is consistent with the phase information present, reconstructs both coding sequences jointly, translates them, and groups arrangements by outcome. It reports whether every arrangement keeps a copy with the reference protein, whether none does, or whether phase decides. It never assigns probabilities, never treats an absent record as reference evidence, and rejects inputs outside its model rather than approximating them.

# Statement of need

Questions such as "are these two variants in *cis* or *trans*?" arise in rare-disease research, carrier analysis and the study of highly polymorphic genes. The answer is often not in the VCF: short-read callers phase only nearby variants, and phase sets from long-read phasers such as WhatsHap [@Martin2016] or LongPhase [@Lin2022] can still end inside a gene. Existing haplotype-aware tools compute consequences for a given phasing; PhaseSpace instead shows the full set of possibilities and what would distinguish them.

PhaseSpace works from local files only: a VCF/gVCF or BCF [@Danecek2011], an indexed reference FASTA and a GTF/GFF3 annotation, read through pysam and HTSlib [@Bonfield2021]. Its outputs are self-contained HTML reports that make no network requests, plus JSON and TSV for downstream use.

# Functionality

- **Single transcript or whole annotation.** `phasespace report` analyses one transcript, chosen by ID or by gene (MANE Select [@Morales2022], then Ensembl canonical, then RefSeq Select). `phasespace scan` analyses one representative transcript per protein-coding gene and ranks genes by whether phase decides if a copy keeps the reference protein.
- **Evidence-aware inputs.** gVCF reference blocks provide explicit coverage of coding bases, subject to genotype-quality and minimum-depth thresholds. Gaps stay visibly unknown. Phase can come from the input or from a separate phased VCF (`--phase-from`).
- **Protein-level description.** Each copy's changes are given in HGVS-style protein notation from its jointly edited sequence. Two calls in one codon therefore give the amino acid actually made: on HG002, two SIRPA calls that alone predict p.Asn100Asp and p.Asn100Lys lie on one copy and give p.Asn100Glu.
- **Read-backed evidence.** With `--reads`, PhaseSpace types each read or read pair at the heterozygous calls it covers and counts fragments supporting the same or opposite orientation for each pair of phase blocks. With `--parents`, the parents' genotypes show which parent passed on each allele. Both kinds of evidence are presented as observations, separate from user hypotheses. They can also link blocks before enumeration (`--link-evidence`), which makes loci with many unphased calls analysable.
- **Interactive exploration.** Reports show both gene copies, protein lengths on a shared scale, per-block flip controls and hypothetical block relationships. They also flag splice-region variants, which are outside the coding model.
- **Workflow integration.** A Docker image and nf-core-style Nextflow modules [@DiTommaso2017; @Ewels2020] are provided.

# Validation

PhaseSpace was validated on the Genome in a Bottle HG002 sample [@Zook2016; @Wagner2022], using GRCh38 chromosome 20.

- **Inputs:** DeepVariant [@Poplin2018] gVCFs from PacBio HiFi and a second DeepTrio run, a GATK HaplotypeCaller [@Poplin2017] gVCF, and the GIAB benchmark VCF.
- **Robustness:** all 521 MANE Select transcripts loaded from every input without error.
- **Consequences:** for 522 benchmark alleles analysed singly, consequence classes agreed with BCFtools/csq in all but one case, a frameshift reading through the stop codon, which PhaseSpace additionally labels as stop loss.
- **Trio phasing:** HG002 was phased by transmission from its parents with WhatsHap pedigree phasing [@Garg2016]. In every transcript, phase only removed arrangements and never created a new outcome.
- **Read evidence:** on HiFi and Illumina reads, orientations agreed with the trio phase wherever both were available, and reads agreed with the trio phase in all 6,969 within-block observations across SIRPA's 19 phased variants.

A 4.7 MB fixture of these real data is part of the test suite.

# Acknowledgements

TODO.

# References
