
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process SNPPILEUP {
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
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path common_vcf

    output:
    tuple val(meta), path("${meta.id}.snppileup.csv.gz"), emit: snppileup_result

    script:
    def prefix     = "${meta.id}.snppileup"

    """
    snp-pileup \\
    --gzip \\
    --min-map-quality 15 \\
    --min-base-quality 20 \\
    --pseudo-snps 100 \\
    --min-read-counts 35 \\
     ${common_vcf} ${prefix}.csv.gz ${normal_bam} ${tumor_bam}
    """
}

process FACETS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/R4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(snppileup_counts)
    path arm_position_rda

    output:
    tuple val(meta), path("${meta.id}.def_cval*.tsv"), path("${meta.id}.def_cval*.txt"), path("${meta.id}.def_cval*.pdf")
    //tuple val(meta), path("${meta.id}.cval*_stats.tsv"), path("${meta.id}.cval*_arm_fga.tsv"), path("${meta.id}.cval*_CNV.txt"), path("${meta.id}.cval*_CNV.pdf"), path("${meta.id}.cval*_CNV_spider.pdf")

    script:
    def snp_nbhd       = 250
    def cval_preproc   = 25
    def cval_proc1     = 25
    //def cval_proc1     = 50
    //def cval_proc2     = 150
    //def cval_proc2     = 100
    def cval_proc2     =50
    def min_read_count = 35
    def ref            = "hg19"
    def mcval          = "NO"
    //def mcval          = "MCVAL"
    def plot           = "PDF"
    """
    Rscript /data/person/huz/pipeline/nf-wgs/bin/facets.cval.r \\
            ${snppileup_counts} \\
            ${ref} ${snp_nbhd} \\
            ${cval_preproc} ${cval_proc1} ${cval_proc2} ${min_read_count}\\
            ${mcval} ${plot}
    """
}

process CNV_FACETS {
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
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_index
    path fasta_dict
    path dbsnp_vcf
    path dbsnp_tbi
    path intervals

    output:
    tuple val(meta), path("b37_cval_1000_500_facets_results/*"), emit: facets

    script:
    def prefix     = "${meta.id}.facets"

    """
    cnv_facets.R \\
    -t ${tumor_bam} \\
    -n ${normal_bam} \\
    -vcf ${dbsnp_vcf} \\
    --cval 500 1000 \\
    --gbuild hg19 \\
    -o b37_cval_1000_500_facets_results/${prefix}
    """
}
