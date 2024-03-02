
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process DOABSOLUTE {
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
    tuple val(meta), path(segment), path(maf)

    output:
    tuple val(meta), path("DoAbsolute/${meta.id}.*"), emit: absolute
    tuple val(meta), path("DoAbsolute/maf/*"), emit: maf
    tuple val(meta), path("DoAbsolute/seg/*"), emit: seg
    tuple val(meta), path("DoAbsolute/summary/*"), emit: summary
    tuple val(meta), path("DoAbsolute/sample_before_summary/*"), emit: sample_before_summary
    tuple val(meta), path("DoAbsolute/sample_final_called/*"), emit: sample_final_called

    script:
    def prefix     = "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(DoAbsolute)
    library(data.table)
    options (warn = -1)
    seg = fread("${segment}")
    colnames(seg) <- c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean")
    seg <- as.data.frame(seg)
    seg [seg[,2] != "Y",]
    maf = fread("${maf}")

    DoAbsolute(Seg = seg, Maf = maf, platform = "Illumina_WES",
           copy.num.type = "total",
           primary.disease = "UCEC",
           results.dir = "DoAbsolute",
           nThread = 2,
           min.mut.af =0.1,
           max.non.clonal =1,
           keepAllResult = TRUE,
           clean.temp = FALSE,
           verbose = FALSE)

   file.rename("DoAbsolute/DoAbsolute.called.ABSOLUTE.plots.pdf","DoAbsolute/${prefix}.called.ABSOLUTE.plots.pdf")
   file.rename("DoAbsolute/DoAbsolute.wsx.ABSOLUTE.table.txt","DoAbsolute/${prefix}.wsx.ABSOLUTE.table.txt")
   file.rename("DoAbsolute/summary/DoAbsolute.PP-modes.plots.pdf","DoAbsolute/summary/${prefix}.PP-modes.plots.pdf")
   file.rename("DoAbsolute/summary/DoAbsolute.PP-calls_tab.txt","DoAbsolute/summary/${prefix}.PP-calls_tab.txt")
   file.rename("DoAbsolute/summary/DoAbsolute.PP-modes.data.RData","DoAbsolute/summary/${prefix}.PP-modes.data.RData")

    """
}
