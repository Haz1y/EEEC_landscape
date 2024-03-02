
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process ASCETS {
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
    tuple val(meta), path(seg)
    path genomic_arm_coor

    output:
    tuple val(meta), path("${meta.id}_aneuploidy_scores.txt"), path("${meta.id}_arm_level_calls.txt"), path("${meta.id}_arm_weighted_average_segmeans.txt"), path("${meta.id}_params.txt"), path("${meta.id}_segmean_hist.pdf")

    script:

    """
    less ${seg} | grep -v 21154917 > ${meta.id}.removeY.seg
    Rscript /data/person/huz/my_software/ascets/run_ascets.R \\
    -i ${meta.id}.removeY.seg -c ${genomic_arm_coor} \\
    -o ${meta.id}
    """
}

process CNVKIT_ASCETS {
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
    tuple val(meta), path(seg)
    path genomic_arm_coor

    output:
    tuple val(meta), path("${meta.id}_aneuploidy_scores.txt"), path("${meta.id}_arm_level_calls.txt"), path("${meta.id}_arm_weighted_average_segmeans.txt"), path("${meta.id}_params.txt"), path("${meta.id}_segmean_hist.pdf")

    script:


    """
    Rscript /data/person/huz/my_software/ascets/run_ascets.R \\
    -i ${seg} -c ${genomic_arm_coor} \\
    -o ${meta.id}
    """
}
