
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process CNVKIT {
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
    path intervals
    path refFlat
    val exome

    output:
    tuple val(meta), path("results/*"), emit: cns

    script:
    def prefix     = "${meta.id}.cnvkit"
    def mode = exome ? "hybrid": " wgs"

    """
    cnvkit.py access \\
    ${fasta} \\
    -o access-excludes.hg19.bed

    cnvkit.py batch \\
    ${tumor_bam} \\
    --normal ${normal_bam} \\
    -m ${mode} \\
    --targets ${intervals} \\
    --annotate ${refFlat} \\
    --fasta ${fasta} \\
    --access ${intervals} \\
    --output-reference ${prefix}.cnn \\
    --output-dir results/ \\
    --diagram --scatter

    """
}
