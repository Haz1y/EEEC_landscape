
conda_env = (params.conda_env ? "${params.conda_env}": "bioconda::fastp")
singularity_url = ""
docker_url = ""

process FASTP {
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
    val   extra_args

    output:
    tuple val(meta), path("*.fastq.gz"),   emit: reads
    tuple val(meta), path("*.fastp.html"), emit: html
    tuple val(meta), path("*.fastp.json"), emit: json

    script:
    def prefix   = "${meta.id}.clean"

    // platform (ILLUMINA DNBSEQ)
    if (meta.platform == "ILLUMINA") {
        adapter_command = "--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    }
    else if (meta.platform == "DNBSEQ") {
        adapter_command = "--adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"
    } else {
        adapter_command = ""
    }

    """
        fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${prefix}.R1.fastq.gz \\
        -O ${prefix}.R2.fastq.gz \\
        ${adapter_command} \\
        ${extra_args} \\
        --thread ${task.cpus} \\
        --html ${prefix}.fastp.html \\
        --json ${prefix}.fastp.json
    """
}
