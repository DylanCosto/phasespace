process PHASESPACE_SCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::phasespace-vcf=0.2.0"
    container "${params.phasespace_container ?: 'ghcr.io/dylancosto/phasespace:0.2.0'}"

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(phased_vcf), path(reads), path(reads_index)
    tuple val(meta2), path(fasta), path(fasta_index)
    path annotation
    path genes

    output:
    tuple val(meta), path("${prefix}/scan.tsv"), emit: tsv
    tuple val(meta), path("${prefix}/scan.json"), emit: json
    tuple val(meta), path("${prefix}/index.html"), emit: html
    tuple val(meta), path("${prefix}/reports"), emit: reports
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample = meta.sample ? "--sample ${meta.sample}" : ''
    def phase = phased_vcf ? "--phase-from ${phased_vcf}" : ''
    def bams = reads ? (reads instanceof List ? reads : [reads]).collect { bam -> "--reads ${bam}" }.join(' ') : ''
    def panel = genes ? "--genes @${genes}" : ''
    """
    phasespace scan \\
        --vcf ${vcf} \\
        --fasta ${fasta} \\
        --annotation ${annotation} \\
        ${sample} ${phase} ${bams} ${panel} \\
        --threads ${task.cpus} \\
        --out-dir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phasespace: \$(phasespace --version | sed 's/PhaseSpace //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/reports
    touch ${prefix}/scan.tsv ${prefix}/scan.json ${prefix}/index.html
    printf '"%s":\\n    phasespace: 0.2.0\\n' "${task.process}" > versions.yml
    """
}
