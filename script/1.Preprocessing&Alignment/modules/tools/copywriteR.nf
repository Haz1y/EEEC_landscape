
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process COPYWRITER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path gc_ref



    output:
    tuple val(meta), path("seg/*"), emit: copywrite_seg

    script:
    def prefix     = "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(CopywriteR)

    bp.param <- SnowParam(workers = ${task.cpus}, type = "SOCK")

    sample.control <-data.frame(samples=c("${tumor_bam}","${normal_bam}"),
                                controls=c("${normal_bam}","${normal_bam}"))

    dir.create("./seg")
    CopywriteR(sample.control,
               destination.folder="./seg",
               reference.folder="./hg19_20kb_chr",
               bp.param,
               keep.intermediary.files = FALSE)

    plotCNA(destination.folder = "./seg")
    unlink("seg/CNAprofiles/BamBaiPeaksFiles",recursive = T, force = FALSE)
    """
}
