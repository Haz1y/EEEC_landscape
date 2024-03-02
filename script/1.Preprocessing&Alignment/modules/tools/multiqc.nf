conda_env = (params.conda_env ? params.conda_env : "bioconda::multiqc=1.10.1")
singularity_url = ""
docker_url = ""

process MULTIQC {
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
    path multiqc_config
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')
    path ('star/*')
    path ('samtools/idxstats/*')
    path ('rseqc/bam_stat/*')
    path ('rseqc/infer_experiment/*')
    path ('rseqc/inner_distance/*')
    path ('rseqc/junction_annotation/*')
    path ('rseqc/junction_saturation/*')
    path ('rseqc/read_distribution/*')
    path ('rseqc/read_duplication/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    script:
    """
    multiqc -f .
    """
}
