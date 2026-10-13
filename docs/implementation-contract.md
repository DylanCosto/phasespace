# Internal contract

Python package `phasespace`, src layout. Coordinates are 0-based half-open.

`inputs.load_locus(vcf_path, fasta_path, annotation_path, transcript_id, sample=None, min_gq=20, min_dp=5) -> Locus`.
Annotation GTF/GFF3. All transcript.cds_intervals in genomic ascending order, reference_cds in transcript orientation. Variant cds_* fields in transcript orientation; genomic ref/alleles in forward orientation. VCF genotype alleles preserved in order. Missing/low-quality calls are evidence unknown, never infer reference from absent records. Supporting optional indexed gVCF blocks. Filter record/sample calls conservatively. Raises PhaseSpaceError with helpful errors for unsupported coding alleles (symbolic called variants, spanning/overlapping complex edits) rather than silently skipping. Plain VCF can stream, indexed files fetch locus.

`analysis.analyze(locus, max_blocks=8) -> dict`. JSON-serializable dict schema:
- schema_version: 1; sample; transcript: asdict(Transcript); variants: asdicts plus block id or null; evidence: asdicts; warnings; source
- blocks: [{id,label,variant_ids: [id]}]. Heterozygotes with same PS share orientation; no PS phased records share implicit block (document); unphased calls independent. Homozygotes fixed on both copies. Preserve 1|0 orientations.
- scenarios: [{id,orientations: {block_id: 0|1},haplotypes: [{variant_ids,cds,protein,consequences:[str],first_stop:int|null,protein_length:int}],outcome_key:str,unaffected_copies:int}]. unaffected_copies means coding sequence identical to supplied reference, not functional. Fix block0 orientation=0 to remove global copy swaps. Never truncate scenarios: exceeding cap errors. Haplotype protein includes translation up to first stop only; first_stop one-based amino-acid position. Include combined edits/frame restoration. Distinct outcomes dedup by unordered pair of full CDS sequences, while retaining scenarios.
- outcome_count; scenario_count; reference_protein; reference_first_stop; reference_protein_length
- informative_pairs: [{block_a,block_b,score:number,description:str}] optional ranking of block orientation parity that separates coding outcomes. score simple count of outcome-pairs distinguished (no probabilistic inference).

`report.write_report(result: dict, output_path) -> None`: standalone interactive HTML with all assets inlined, no remote requests. Render phase constraints from orientations. User constraints are explicitly hypothetical; never mutate input evidence. Show sample/transcript, evidence strips, variants, two copy tracks, predicted proteins, coding consequence labels, result counts, reset constraints, JSON export. HTML escape via safe JSON embedding and textContent. Be responsive/accessibility conscious. Do not claim biological function or diagnosis. Plain VCF absence is unknown.

Root owns models.py, CLI, packaging, docs, synthetic examples, integration tests. Separate agents own inputs.py + tests/test_inputs.py; analysis.py + tests/test_analysis.py; report.py/assets + tests/test_report.py. No dependency manifest edits by agents.
