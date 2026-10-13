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

## Important scope

These tests establish behavior for the supported computational model. They do not establish biological function or clinical performance. Transcript annotation completeness, genotype quality, the correct reference build and supported genetic code remain necessary inputs. Missing evidence never becomes an inferred reference call.
