
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process PYCLONE {
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
    tuple val(meta), path(c_input), path(ln_input)

    output:
    tuple val(meta), path("output/*"), emit: output
    tuple val(meta), path("log/*"), emit: log

    script:
    def prefix     = "${meta.id}"

    """
    mkdir output
    mkdir log

    PyClone run_analysis_pipeline \\
    --prior total_copy_number \\
    --in_files ${c_input} ${ln_input} \\
    --working_dir ./output/${prefix}_pyclone_analysis 1>./log/${prefix}.log 2>&1
    """
}
