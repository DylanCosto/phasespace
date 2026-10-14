# PhaseSpace for Nextflow

`modules/phasespace/scan` and `modules/phasespace/report` are nf-core-style modules; `main.nf` is a small example pipeline that runs one scan per sample.

```bash
nextflow run nextflow/main.nf -profile docker \
  --samplesheet samples.csv --fasta GRCh38.fa --annotation gencode.v46.annotation.gtf.gz \
  [--genes panel.txt] [--outdir results]
```

`samples.csv` columns: `sample`, `vcf`, and optionally `vcf_index` (default `<vcf>.tbi`), `vcf_sample` (for multi-sample VCFs), `phased_vcf` (passed as `--phase-from`) and `reads` (BAM/CRAM paths separated by `;`, each with an index). Relative paths are resolved against the samplesheet's folder. The reference needs a `.fai` (and `.gzi` if bgzipped).

Profiles: `docker`, `singularity` and `conda`. Set `--phasespace_container` to use another image; the default is `ghcr.io/dylancosto/phasespace:<version>`. Add extra arguments per process with `ext.args`, e.g. `withName: PHASESPACE_SCAN { ext.args = '--reports all --max-blocks 10' }`.

Try it on the bundled HG002 fixture:

```bash
nextflow run nextflow/main.nf -profile docker,test
```

To use a module in your own pipeline, copy `modules/phasespace/<name>` and include it:

```groovy
include { PHASESPACE_SCAN } from './modules/phasespace/scan/main'
PHASESPACE_SCAN(samples, reference, annotation, genes)   // pass [] for optional inputs
```
