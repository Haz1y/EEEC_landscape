
conda_env = (params.conda_env? params.conda_env : "bioconda::bwa=0.7.17 bioconda::samtools=1.12")
singularity_url = ""
docker_url = ""

process BWA_MEM {
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
    path  index
    val   extra_args

    output:
    tuple val(meta), path("*.mapped_sorted.bam"), emit: bam

    script:
    def cpus = task.cpus
    def prefix     = "${meta.id}"
    def read_group = "-R \"@RG\\tID:${prefix}\\tSM:${meta.sample}\\tPL:${meta.platform}\""
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
        $extra_args \\
        $read_group \\
        -t ${cpus} \\
        \$INDEX \\
        $reads \\
        | samtools view -@ ${cpus} -bhS -o ${prefix}.mapped.bam -
    samtools sort --threads ${cpus} -m 2G ${prefix}.mapped.bam -o ${prefix}.mapped_sorted.bam
    rm -rf ${prefix}.mapped.bam
    """
}
