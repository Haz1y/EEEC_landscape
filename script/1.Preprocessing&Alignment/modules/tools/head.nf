
conda_env = (params.conda_env? params.conda_env : "bioconda::bwa=0.7.17 bioconda::samtools=1.12")
singularity_url = ""
docker_url = ""

process HEAD {
    tag "$meta.id"
    label 'process_high'
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
    tuple val(meta), path("*.r{1,2}.fastq"), emit: fastq

    script:
    def prefix     = "${meta.id}"
    """
    zcat ${reads[0]} | head -n 40000 > ${prefix}.r1.fastq
    zcat ${reads[1]} | head -n 40000 > ${prefix}.r2.fastq
    """
}
