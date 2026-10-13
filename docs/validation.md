# Validation record — version 0.1.0

Validated locally on 2026-10-12 (America/New_York), Linux x86_64, Python 3.12.14, pysam 0.24.1, and Node 26.7.0.

- **107 automated tests passed.** This includes 55 input tests, 42 analysis tests, six command-line integration tests, and four report tests. Report tests also execute JavaScript filtering checks with Node.
- Browser checks exercised unresolved, same-orientation and opposite-orientation hypotheses. The two-stop synthetic example produced stop positions 8/32 in cis and 8/21 in trans, matching the computed JSON.
- The example gVCF displayed 88 reference-called bases, two variant-called bases and six unknown bases within its 96-base CDS. The corresponding variant-only VCF left 94 bases unknown.
- The browser JSON download was read back successfully, with the original analysis separated from hypothetical exploration state.
- Desktop and 390-pixel viewport layouts were inspected. The narrow layout had no document-level horizontal overflow. Browser error logs were empty during these checks.
- Source and wheel packages built successfully. The installed wheel generated its demo outside the source checkout, checking that inline report assets were packaged.

## Speed check

`python scripts/benchmark.py` creates eight independent heterozygous calls in a synthetic 3,000-base CDS and evaluates all 128 arrangements. Five local runs had a median engine time of **0.3115 seconds**, with a range of **0.2446–0.3639 seconds**. This excludes input parsing and report generation. It is not a whole-genome or clinical dataset benchmark.

The synthetic demo's entire locus analysis, including indexed input loading, took approximately 0.008 seconds in a warmed local run. Real annotation scans and file access will vary substantially. Large real datasets, diverse callers, and additional organisms have not yet been benchmarked or independently validated.

## Update — 2026-10-13

Checked on Linux x86_64 with Python 3.12, pysam 0.24.1 and Node 26.

- **117 automated tests passed**, including new checks for gVCF deletions (reference blocks through deleted bases, `*` spanning-deletion alleles, low-quality spanning records, and rejection of `*` combined with another alternate allele). Report tests check the compact encoding, confirm in Node that the browser rebuilds the exact analysis from it, and cover the phase verdict, block flipping and coding-change labels.
- `python scripts/benchmark.py` measured a median engine time of **0.043 s** on this machine, down from 0.144 s before the translation changes on the same machine. The figures in the section above came from a different machine.
- At the 12-block cap with a 10,002-base CDS (2,048 arrangements), analysis takes about 3 s, down from 9.6 s. The HTML report is **1.6 MB**, down from 44 MB, and loads in about 130 ms in a desktop browser.
- The redesigned report was inspected at desktop and 375-pixel widths, in light and dark modes, with no document-level horizontal overflow and no console errors.

## Real data — GIAB HG002, GRCh38 chr20 (2026-10-13)

Inputs (all public): DeepVariant 1.10 DeepTrio gVCFs for HG002 (a PacBio chr20 run, and a second run published under `ont-r104-chr20` whose file also holds exome-scale calls elsewhere, so its platform is not asserted), the GIAB NISTv4.2.1 benchmark VCF, GENCODE v46 and Ensembl 112 annotation, and chr20 of the GRCh38 no-alt analysis set.

- **Sweep.** Every chr20 MANE Select transcript (521) was run against four inputs (both gVCFs, the DeepVariant variant-only VCF and GIAB): 2,084 runs. All loaded without error. Four stopped at the 12-block cap with 20–21 heterozygous coding variants, as designed. The callers agreed on the variants in the transcripts inspected closely.
- **Annotation.** With GTF `stop_codon` features joined, all 521 GENCODE coding sequences end in a stop codon and match the Ensembl GFF3 intervals exactly, whether the transcript is named with or without its version. Before this change, every GENCODE run warned that the terminal stop was missing and could not assess stop loss.
- **Consequences.** For 522 GIAB variant alleles in MANE coding sequence, each analysed alone, PhaseSpace agreed with `bcftools csq --local-csq` (Ensembl 112) on every consequence class (missense, synonymous, stop gained or lost, frameshift, in-frame insertion or deletion, start lost) except one. That exception is the common DEFB126 frameshift that reads through the stop codon: PhaseSpace reports `frameshift` and `stop_lost`, bcftools reports only `frameshift`.
- **Phase verdict.** 49 transcripts with both gVCFs had a reference-identical CDS in only some arrangements. In 19 of these, every change was synonymous, so the report verdict now asks whether a copy keeps the reference protein sequence.
- **gVCF deletions.** DeepVariant 1.10 starts the next reference block after a deletion's REF span, so these files never needed the overlapping-block handling. That handling remains for writers that do overlap, and is covered by synthetic tests.
- **Speed.** Looking up one transcript in the full 51 MB GENCODE v46 GTF took 4.2 s, down from 59 s, by skipping lines that cannot belong to it before parsing them.
- **Trio phasing.** None of the public files above carries phased genotypes, so HG002 was phased by transmission. Its DeepTrio gVCF was merged with the HG003 and HG004 gVCFs from the same run (`bcftools merge --gvcf`), genotypes with GQ < 20 were set missing, and the result was phased with WhatsHap 2.8 (`phase --ped`, no reads). Of HG002's heterozygous chr20 calls, 1,637 were phased into one chromosome-wide phase set and 882 stayed unphased (mostly sites heterozygous in all three samples). `scripts/transfer_phase.py` copied the phase back onto the gVCF, matching 1,636 records by allele sequence. Phase was not checked against an independent phasing truth set.
- **Phased versus unphased.** Running all chr20 MANE transcripts on the unphased and phased versions of the same gVCF, every phased outcome was among the unphased outcomes (520 transcripts, no exceptions). Phase reduced the arrangements in 61 transcripts. Of 17 transcripts where phase decided whether a copy keeps the reference protein, phasing resolved 14: in 8 a copy keeps it, in 6 neither does. The other 3 retain unphased sites. SIRPA, with 20 heterozygous coding variants, exceeds the 12-block cap unphased but has 2 arrangements when phased: one copy is reference-identical in its CDS, and the other carries 19 alternate alleles, including an in-frame deletion.
- **Report on real phased data.** SIRPA's 19 clustered variants exposed overlapping pin labels. Labels are now limited to three rows, never cross a stem, and any that do not fit are listed below the track.
- **Regression fixture.** `tests/data/hg002_chr20` (424 KB, made by `scripts/make_hg002_fixture.py`) keeps real records for six of these transcripts (five plus SIRPA), including the trio-phased gVCF, in true chr20 coordinates. `tests/test_hg002.py` checks each across three inputs and both annotation formats.

## Important scope

These tests establish behavior for the supported computational model. They do not establish biological function or clinical performance. Transcript annotation completeness, genotype quality, the correct reference build and supported genetic code remain necessary inputs. Missing evidence never becomes an inferred reference call.
