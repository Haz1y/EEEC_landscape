
conda_env = (params.conda_env? params.conda_env : "bioconda::samtools=1.12")
singularity_url = ""
docker_url = ""

process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.merged.bam"), emit: bam

    script:
    prefix   = "${meta.id}"
    """
    samtools merge ${prefix}.merged.bam ${bams}
    """
}
