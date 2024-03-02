
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process PURBAYES {
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
    tuple val(meta), path(maf)


    output:
    tuple val(meta), path("${meta.id}.PurBayes.*.csv"), emit: purbayes_out
    tuple val(meta), path("${meta.id}.PurBayes.plot.pdf"), emit: purbayes_plot

    script:
    def prefix     = "${meta.id}.PurBayes"

    """
    #!/usr/bin/env Rscript
    library(readr)
    library(PurBayes)

    maf = read_tsv("${maf}",skip=1)
    maf = as.data.frame(maf)
    t_depth = maf[,"t_depth"]
    t_alt = maf[,"t_alt_count"]
    Pur = PurBayes(t_depth,t_alt)
    out = summary(Pur)

    write.table(t(out[[1]]),sep=",",file="${prefix}.purity.csv",quote=F,row.names=F)
    write.table(cbind(out[[2]][,1],out[[2]]),sep=",",file="${prefix}.post.dist.csv",quote=F,row.names=T)
    write.table(cbind("${meta.id}",out[[3]]),sep=",",file="${prefix}.n.pop.csv",quote=F,col.names=F,row.names=F)
    write.table(out[[4]],sep=",",file="${prefix}.dev.out.csv",quote=F,row.names=F)

    pdf("${prefix}.plot.pdf")
    plot(Pur)
    dev.off()

    """
}
