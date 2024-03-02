
process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    input:
    tuple val(meta), path(bam), path(bai)
    val intervals

    output:
    tuple val(meta), path("${meta.id}_bamqc_result"), emit: bamqc_report

    script:
    def prefix   = "${meta.id}.manta"
    def intervals_cmd  = intervals ? "--exome --callRegions ${intervals}" : ""
    """

    /data/person/huz/my_software/qualimap-build-31-08-20/qualimap bamqc \\
    -bam ${bam} \\
    --java-mem-size=${task.memory.toGiga()}G \\
    -nt ${task.cpus} \\
    -outdir ${meta.id}_bamqc_result \\
    -outfile ${meta.id}_bamqc_report.pdf \\
    -outformat PDF:HTML \\
    -gff ${intervals}

    """
}
