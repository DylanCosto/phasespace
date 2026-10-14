# Model and interpretation

## Coordinates and transcript sequence

Internal genomic and CDS coordinates are zero-based, half-open. VCF and annotation display positions are one-based. CDS segments are stored in ascending genomic order; the joined reference sequence is reverse-complemented for a negative-strand transcript. Internal CDS phase bases are retained because they finish codons split across exon boundaries. Input phase values validate the joining rather than trim every segment.

Every overlapping VCF REF allele is checked against the supplied FASTA. A sequence-resolved variant must fit within one CDS segment before it is projected into transcript coordinates. This deliberately bounds the model instead of guessing a splicing consequence.

## Phasing

For each heterozygous call, the two allele sequences retain the order encoded in its genotype. Calls with a common `PS` value and phased separator belong to one orientation block. Phased calls without `PS` share an implicit block. Each unphased heterozygote is independent; `PS` on an unphased genotype does not link it to other calls. Homozygous calls are applied to both copies in every arrangement.

The first block's orientation is fixed to remove the symmetry from renaming the two copies. With B independent blocks this leaves 2^(B−1) arrangements when B is positive, and one when B is zero. All supported arrangements are evaluated. The number of blocks is bounded before enumeration.

Edits are applied jointly to the coding sequence. Translation uses the standard genetic code and ends at the first stop. Different arrangements can yield the same unordered pair of complete coding sequences; they remain inspectable while sharing an outcome key. Distinct coding outcomes can still have identical protein sequences. These labels are not statistical confidence intervals or posterior probabilities.

## Evidence

Reference evidence requires an explicit usable homozygous-reference call or reference block. A site absent from a variant-only VCF provides no reference evidence. Failed calls and missing quality fields remain unknown. Quality thresholds apply to the selected sample. A multi-base block needs `MIN_DP` for a positive depth threshold; its anchor `DP` is not used as evidence of minimum depth throughout the interval. A zero threshold explicitly waives that block-wide depth check.

DeepVariant writes confident homozygous-reference calls at candidate sites with `FILTER=RefCall`. On a `0/0` genotype that filter is treated as passing, so the call counts as reference evidence if it meets the GQ and depth thresholds. `RefCall` on any other genotype still excludes the call.

A called deletion carries variant evidence across its entire REF span. Some gVCF writers continue reference blocks through the bases a deletion removes from one copy, so reference calls over those removed bases are dropped instead of being reported as conflicting. The deletion's anchor base gets no such exemption: a passing reference call there remains a conflict. A `*` spanning-deletion record restates an upstream deletion. It adds no edit, and its own quality still decides whether its span is variant or unknown evidence.

Evidence is partitioned over CDS intervals. Introns are not counted as missing coding evidence. Unknown evidence does not manufacture additional variants; the reported arrangements are conditional on the included calls. All sequence-level interpretations must be read with that qualification.

## Hypothesis controls

The report's controls constrain the **orientation of entire phase blocks**. They are explicitly hypothetical. Two variants already oriented oppositely within one block stay opposite when the block flips. Consequently, aligning two blocks does not necessarily put every alternate allele from those blocks on the same copy.

Pair rankings are an exploration aid for distinguishing model outcomes. They are not a recommendation for a particular laboratory assay, a proof of functional effect, or a probabilistic estimate of which arrangement is true.

## Validation strategy

Tests cover small cases with known outcomes, including joint codon changes, premature stops in cis/trans, mixed genotype orientations, homozygous and multi-allelic calls, frame-restoring edits, reverse-strand projection, reference evidence gaps, quality failures, and unsupported-input rejection. End-to-end tests use a generated synthetic FASTA/GTF/VCF/gVCF. Release packaging must include all report assets and the installed demo must work outside the checkout.
