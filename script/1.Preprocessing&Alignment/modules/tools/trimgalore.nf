conda_env = (params.conda_env ? params.conda_env : "bioconda::trim-galore=0.6.6")
singularity_url = ""
docker_url = ""

process TRIMGALORE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode
        enabled: params.publish_enable

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz")    , emit: reads
    tuple val(meta), path("*report.txt"), emit: log
    tuple val(meta), path("*.html"), emit: html optional true
    tuple val(meta), path("*.zip") , emit: zip optional true

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    // --clip_R1 7
    // --clip_R2 7
    def prefix   = "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        trim_galore \\
            --cores $cores \\
            --gzip \\
            ${prefix}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${prefix}.R1.fastq.gz ] && ln -s ${reads[0]} ${prefix}.R1.fastq.gz
        [ ! -f  ${prefix}.R2.fastq.gz ] && ln -s ${reads[1]} ${prefix}.R2.fastq.gz
        trim_galore \\
            --fastqc \\
            --cores $cores \\
            --paired \\
            --gzip \\
            ${prefix}.R1.fastq.gz \\
            ${prefix}.R2.fastq.gz
        """
    }
}
