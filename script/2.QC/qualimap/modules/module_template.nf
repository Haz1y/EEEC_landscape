/*
Template of nextflow process for processing HTS data

*/

conda_env = (params.conda_env? params.conda_env : "")
singularity_url = ""
docker_url = ""

process HELLOWORLD {
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
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    def prefix     = "${meta.id}"

    """
    bwa mem \\
        ${index} \\
        ${reads} \\
        > ${prefix}.sam

    """
}
