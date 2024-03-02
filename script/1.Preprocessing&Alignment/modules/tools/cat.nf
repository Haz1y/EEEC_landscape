
conda_env = (params.conda_env ? params.conda_env : "conda-forge::sed=4.7")
singularity_url = ""
docker_url = ""

process CAT_FASTQ {
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

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def prefix   = "${meta.id}"
    def readList = reads.collect{ it.toString() }
    if (meta.single_end) {
        if (readList.size > 1) {
            """
            cat ${readList.sort().join(' ')} > ${prefix}.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fastq.gz
            """
        }
    }
}
