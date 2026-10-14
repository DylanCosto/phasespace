#!/usr/bin/env nextflow
// Example pipeline: one PhaseSpace scan per sample in a samplesheet.
// nextflow run nextflow/main.nf -profile docker --samplesheet samples.csv --fasta ref.fa --annotation genes.gtf.gz

include { PHASESPACE_SCAN } from './modules/phasespace/scan/main'

params.samplesheet = null
params.fasta = null
params.annotation = null
params.genes = null
params.outdir = 'results'

// Samplesheet paths may be absolute or relative to the samplesheet's folder.
def sheetFile(String value, base) {
    return value.startsWith('/') || value.contains('://') ? file(value, checkIfExists: true) : file(base.resolve(value).toString(), checkIfExists: true)
}

def optionalFile(value) {
    return value ? file(value, checkIfExists: true) : []
}

workflow {
    if (!params.samplesheet || !params.fasta || !params.annotation) {
        error "Provide --samplesheet, --fasta and --annotation"
    }
    def base = file(params.samplesheet, checkIfExists: true).parent
    samples = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def vcf = sheetFile(row.vcf, base)
            def index = row.vcf_index ? sheetFile(row.vcf_index, base) : file("${vcf}.tbi", checkIfExists: true)
            def phased = row.phased_vcf ? sheetFile(row.phased_vcf, base) : []
            def bams = row.reads ? row.reads.tokenize(';').collect { bam -> sheetFile(bam, base) } : []
            def bais = bams.collect { bam -> file("${bam}.bai").exists() ? file("${bam}.bai") : file("${bam}.crai") }
            [[id: row.sample, sample: row.vcf_sample ?: null], vcf, index, phased, bams, bais]
        }
    def fasta = file(params.fasta, checkIfExists: true)
    def fastaIndex = [file("${params.fasta}.fai", checkIfExists: true)] + (file("${params.fasta}.gzi").exists() ? [file("${params.fasta}.gzi")] : [])
    reference = channel.value([[id: fasta.baseName], fasta, fastaIndex])
    PHASESPACE_SCAN(samples, reference, file(params.annotation, checkIfExists: true), optionalFile(params.genes))
}
