
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process SCARHRD {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta), path(small_seqz), path(small_seqz_index)



    output:
    tuple val(meta), path("*_HRDresults.txt"), path("*_info_seg.txt"),emit: hrd_results

    script:

    """
    #!/usr/bin/env Rscript
    library(scarHRD)
    scar_score("${meta.id}.small.seqz.gz",
                chr.in.names =TRUE,
	             reference = "grch37", seqz=FALSE)
    """
}
