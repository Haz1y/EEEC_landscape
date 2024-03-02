
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process MSISENSOR {
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
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_index
    path fasta_dict
    path ref_site

    output:
    tuple val(meta), path("results/*"), emit: cns

    script:
    def prefix     = "${meta.id}.msisensor"

    """
    msisensor-pro msi \\
    -d ${ref_site} \\
    -n ${normal_bam} \\
    -t ${tumor_bam} \\
    -o results/${prefix}
    """
}
