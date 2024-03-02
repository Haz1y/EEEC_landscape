conda_env = (params.conda_env ? "${params.conda_env}":"bioconda::fastqc=0.11.9")
singularity_url = ""
docker_url = ""

process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
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
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix   =  "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fq.gz ] && ln -s $reads ${prefix}.fq.gz
        fastqc --threads $task.cpus ${prefix}.fq.gz
        """
    } else {
        """
        [ ! -f  ${prefix}.R1.fastq.gz ] && ln -s ${reads[0]} ${prefix}.R1.fastq.gz
        [ ! -f  ${prefix}.R2.fastq.gz ] && ln -s ${reads[1]} ${prefix}.R2.fastq.gz
        fastqc --threads $task.cpus ${prefix}.R1.fastq.gz ${prefix}.R2.fastq.gz
        """
    }
}
