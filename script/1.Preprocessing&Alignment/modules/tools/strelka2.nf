
conda_env = (params.conda_env? params.conda_env : "bioconda::strelka=2.9.10")
singularity_url = "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
docker_url = "quay.io/biocontainers/strelka:2.9.10--0"

process STRELKA2_GERMLINE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable


    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq_py2.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    val intervals
    val extra_args

    output:
    tuple val(meta), path("*.variants.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*.genome.vcf.gz")      , emit: genome_vcf
    tuple val(meta), path("*.genome.vcf.gz.tbi")  , emit: genome_vcf_tbi

    script:
    def prefix   = "${meta.id}"
    def intervals_cmd  = intervals ? "--exome --callRegions ${intervals}" : ""
    """
    configureStrelkaGermlineWorkflow.py \\
        --bam ${bam} \\
        --referenceFasta ${fasta} \\
        ${intervals_cmd} \\
        ${extra_args} \\
        --runDir ./

    python ./runWorkflow.py -m local -j $task.cpus
    mv ./results/variants/genome.*.vcf.gz     ${prefix}.genome.vcf.gz
    mv ./results/variants/genome.*.vcf.gz.tbi ${prefix}.genome.vcf.gz.tbi
    mv ./results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
    mv ./results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi
    """
}


process STRELKA2_SOMATIC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq_py2.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(candi_indel_vcf), path(candi_indel_tbi)
    path  fasta
    path  fai
    path  dict
    path intervals

    output:
    tuple val(meta), path("*.somatic.snvs.vcf.gz"), path("*.somatic.snvs.vcf.gz.tbi"), path("*.somatic.indels.vcf.gz"), path("*.somatic.indels.vcf.gz.tbi") , emit: snv_indels

    script:
    def prefix   = "${meta.id}.strelka"
    def exome_cmd  = params.exome ? "--exome" : ""
    """
    bgzip --threads ${task.cpus} -c ${intervals} > call_targets.bed.gz ; tabix call_targets.bed.gz


    configureStrelkaSomaticWorkflow.py \
        --tumor ${tumor_bam} \
        --normal ${normal_bam} \
        --referenceFasta ${fasta} \
        --indelCandidates ${candi_indel_vcf} \
        --callRegions call_targets.bed.gz \
        ${exome_cmd} \
        --runDir .

    python ./runWorkflow.py -m local -j $task.cpus
    mv ./results/variants/somatic.indels.vcf.gz     ${prefix}.somatic.indels.vcf.gz
    mv ./results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic.indels.vcf.gz.tbi
    mv ./results/variants/somatic.snvs.vcf.gz     ${prefix}.somatic.snvs.vcf.gz
    mv ./results/variants/somatic.snvs.vcf.gz.tbi ${prefix}.somatic.snvs.vcf.gz.tbi
    """
}
