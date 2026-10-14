# Architecture

PhaseSpace is a Python package in `src/phasespace` with three stages and a CLI. Every stage passes plain data forward, so each one can be tested on its own. Internal coordinates are zero-based and half-open. Positions shown to users are one-based.

```
inputs.load_locus  ──▶  Locus  ──▶  analysis.analyze  ──▶  result dict  ──▶  report.write_report  ──▶  HTML
   (VCF/gVCF, FASTA,                (enumerate phase                         (inline CSS/JS, data
    GTF/GFF3, phase VCF)             arrangements)                            embedded as JSON)
```

## `models.py`

Frozen dataclasses shared by every stage: `Transcript`, `Variant`, `Evidence`, `Locus`. `PhaseSpaceError` is the one error type for anything outside the supported model. The CLI prints it and exits with status 2.

## `annotation.py` — finding transcripts

`read_transcript(path, id)` finds one transcript (exact, versioned or unversioned, or GFF3 ID). `select_transcript(path, gene)` picks a gene's MANE Select, then Ensembl canonical, then RefSeq Select transcript. `catalogue(path, genes=None)` returns one such transcript per protein-coding gene in two passes: one over transcript features and one over CDS and stop codons. All three build intervals with the same checks.

## `inputs.py` — reading one locus

`LocusReader(vcf, fasta, sample, …, phase_from)` opens the inputs once and `load(model)` builds each `Locus`. `scan` relies on this.


`load_locus(vcf_path, fasta_path, annotation_path, transcript_id, sample=None, min_gq=20, min_dp=5, phase_from=None) -> Locus`

- **Annotation.** The transcript's `CDS` features define the coding sequence. GTF `stop_codon` features are joined to it, because GTF leaves the stop codon out of `CDS`. A transcript ID resolves exactly; or as a GFF3 `ID`, `transcript_id` or `transcript_id.version`; or without its version when only one annotated version exists. Lines that cannot belong to the transcript are skipped before parsing. Partial models, inconsistent CDS phases and translation exceptions are rejected.
- **Contigs.** The annotation, FASTA, VCF and phase file may name a contig with or without a `chr` prefix (`chrM` and `MT` are also equivalent). User-facing positions and variant IDs use the VCF's name.
- **Evidence.** Reference evidence needs an explicit passing homozygous-reference call or a gVCF block. A missing record is unknown evidence, never reference. DeepVariant's `FILTER=RefCall` on a `0/0` genotype counts as reference, subject to the same GQ and depth thresholds. Multi-base blocks need `MIN_DP` when the depth threshold is positive. A called deletion's REF span is variant evidence, and reference blocks running through its deleted bases are not treated as conflicts.
- **Calls.** Alleles are projected into transcript orientation. Variant IDs are `contig:pos:REF>ALT[,ALT]`, with the alternate alleles sorted so the ID does not depend on genotype order. Symbolic alleles, overlapping edits, edits crossing CDS boundaries and REF mismatches raise `PhaseSpaceError`. A `*` allele is accepted only when a called upstream deletion covers it.
- **Phase.** Genotype order and `PS` come from the input; phased calls without `PS` share an implicit set. With `phase_from`, phase comes only from that VCF. Calls match by position, REF and allele sequence, and any heterozygous call the phase file does not phase becomes unphased.

`check_thresholds(min_gq, min_dp)` validates thresholds; the CLI calls it before reading any input.

## `analysis.py` — enumerating arrangements

`analyze(locus, max_blocks=8) -> dict` (JSON-serialisable, `schema_version` 1)

- Heterozygous calls group into blocks: one per explicit `PS`, one implicit block for phased calls without `PS`, and one per unphased call. Homozygous calls sit on both copies.
- The first block's orientation is fixed, removing the copy-label symmetry, so B blocks give 2^(B−1) arrangements. Exceeding `max_blocks` (0–12; `check_max_blocks`) raises an error. Arrangements are never sampled or truncated.
- Each arrangement (`scenarios[]`) has two copies (`haplotypes[]`). Each copy records `variant_ids`, the `alleles` it carries, `cds`, `protein` (up to the first stop), `consequences`, `first_stop` and `protein_length`. Consequences use the standard genetic code and track base origins, so in-frame indels do not create false stop gains or losses.
- `outcome_key` groups arrangements by their unordered pair of CDS sequences, and `protein_outcome_key` by their predicted proteins. `informative_pairs` ranks block pairs by how many protein outcomes their relative orientation separates. This is not a probability.

## `report.py` and `assets/` — the HTML report

`write_report(result, path)` inlines `report.css` and `report.js` into `report.html` and embeds the result as JSON. Strings are escaped so data can never become markup, and the page makes no network requests (CSP `connect-src 'none'`).

To keep files small, the embedded result drops each copy's `cds` and `protein` (`report_encoding: "edits"`). `report.js` rebuilds them from the reference and the copy's alleles when needed, and the JSON export restores the full schema. The browser-side functions (`filterScenarios`, `expandAnalysis`, `phaseVerdict`, `flipTarget`, `changeLabel`, `placeLabels`, …) are exported for Node, and `tests/test_report.py` tests them directly.

## `scan.py` — whole-annotation scans

`run_scan` takes the catalogue, drops genes on contigs absent from the inputs, and loads and analyses each gene in worker processes. Each worker keeps one `LocusReader` open. `analysis.summarize` gives the protein-level verdict and tier, matching the report's headline; a Node test checks that they agree. Outputs are `scan.tsv`, `scan.json`, `index.html` (from `assets/scan.*`, with the same safety rules as reports) and `reports/`.

## `reads.py` — read-backed evidence

`read_evidence(result, bams, fasta)` types each read at the heterozygous calls it covers: SNVs by aligned base, indels by exact match of a flanked window to either expected allele. Typings join into fragments by read name. For each fragment and pair of calls, allele indices i₁ and i₂ (genotype order) give relative orientation i₁ XOR i₂. Block pairs count fragments by majority vote; same-block pairs check the input phase. `add_read_evidence` adds `supported` relationships (≥ 3 fragments, ≥ 90% agreement) and the verdict among the arrangements they allow. The report shows the counts and can apply them as hypotheses tagged "reads".

## `hgvs.py` — c. and p. labels

`coding_change` gives c. labels (exact for substitutions; indels in VCF representation). `protein_changes` describes one copy from its joint sequence. It walks reference and alternate codons using the alignment from `analysis._alignment`, and names substitutions, nonsense, in-frame del/ins/dup/delins (trimmed to the changed residues), frameshifts with the distance to the new stop, stop loss and start loss. `analyze` stores these per copy (`protein_changes`) and, for each variant on its own, in `variants[].effects`.

## `inheritance.py` — parent-of-origin evidence

For each heterozygous child call, the parents' genotypes (at the site, or from a covering gVCF reference block; filters and GQ required, depth when present) decide which allele came from the father when only one assignment fits. Calls inherited from the same parent share a copy, so two calls give relative orientation f₁ XOR f₂. Evidence has the same shape as read evidence and shares `reads.apply_supported`. A relationship counts as supported when every informative pair agrees.

## `linking.py` — `--link-evidence`

Before enumeration, clearly supported read and inheritance relationships merge blocks through a parity union-find, strongest first. Relationships that contradict stronger ones are skipped and reported. Linked calls are rewritten as one phase set, so loci over the block cap can be analysed, and labels and warnings say the arrangements assume this evidence.

## Tests and fixtures

- `tests/test_inputs.py`, `test_analysis.py`, `test_report.py`, `test_cli.py`: synthetic cases with known answers.
- `tests/test_hg002.py`: real GIAB HG002 chr20 calls in `tests/data/hg002_chr20`, made by `scripts/make_hg002_fixture.py`. They cover six transcripts, three callers or inputs, trio-phased data and both annotation formats. See `docs/validation.md` for how the data were obtained and checked.
