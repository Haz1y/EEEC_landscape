conda_env = (params.conda_env ? params.conda_env : 'bioconda::star=2.6.1d bioconda::samtools=1.10')
singularity_url = ""
docker_url = ""

process STAR_ALIGN {
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
    path  index


    output:
    tuple val(meta), path("*.sorted.bam"),      emit: bam
    tuple val(meta), path("*.sorted.bam.bai"),      emit: bai
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript

    script:
    def prefix     =  "${meta.id}.star"
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --readFilesCommand zcat \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --outSAMtype BAM Unsorted \\
        --quantMode TranscriptomeSAM GeneCounts
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam ${prefix}.*d.out.bam
    samtools index ${prefix}.sorted.bam
    rm -rf ${prefix}.Aligned.out.bam
    """
}
