process PHASESPACE_REPORT {
    tag "${meta.id}:${gene}"
    label 'process_single'

    conda "bioconda::phasespace-vcf=0.2.0"
    container "${params.phasespace_container ?: 'ghcr.io/dylancosto/phasespace:0.2.0'}"

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(phased_vcf), path(reads), path(reads_index), val(gene)
    tuple val(meta2), path(fasta), path(fasta_index)
    path annotation

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val(meta), path("${prefix}.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${gene}"
    def sample = meta.sample ? "--sample ${meta.sample}" : ''
    def phase = phased_vcf ? "--phase-from ${phased_vcf}" : ''
    def bams = reads ? (reads instanceof List ? reads : [reads]).collect { bam -> "--reads ${bam}" }.join(' ') : ''
    """
    phasespace report \\
        --vcf ${vcf} \\
        --fasta ${fasta} \\
        --annotation ${annotation} \\
        --gene ${gene} \\
        ${sample} ${phase} ${bams} \\
        --out ${prefix}.html \\
        --json ${prefix}.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phasespace: \$(phasespace --version | sed 's/PhaseSpace //')
    END_VERSIONS
    """
}
