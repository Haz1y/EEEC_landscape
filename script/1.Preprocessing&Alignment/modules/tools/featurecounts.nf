conda_env = (params.conda_env ? params.conda_env : "bioconda::subread=2.0.1")
singularity_url = ""
docker_url = ""

process SUBREAD_FEATURECOUNTS {
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
    tuple val(meta), path(bams)
    path(gtf)
    path(bed)
    val is_saf

    output:
    tuple val(meta), path("*counts.tsv")        , emit: counts
    tuple val(meta), path("*featureCounts.txt")        , emit: results
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary

    script:
    def prefix     = "${meta.id}"
    def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }

    index = gtf
    saf_command = ""
    if (is_saf) {
      index = bed
      saf_command = "-F SAF"
    }

    """
    featureCounts \\
        $paired_end \\
        -T $task.cpus \\
        -a $index \\
        -g gene_id \\
        -t exon \\
        -s $strandedness \\
        $saf_command \\
        -B \\
        -C \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}

    cut -f 1,7 ${prefix}.featureCounts.txt | tail -n +3 > ${prefix}.counts.tsv
    """
}
