
conda_env = (params.conda_env? params.conda_env : "bioconda::strelka=2.9.10")
singularity_url = "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
docker_url = "quay.io/biocontainers/strelka:2.9.10--0"


process MANTA_GERMLINE {
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
    tuple val(meta), path("*.diploidSV.vcf.gz"), path("*.diploidSV.vcf.gz.tbi") , emit: diploid_vcf
    tuple val(meta), path("*.candidateSV.vcf.gz"), path("*.candidateSV.vcf.gz.tbi") , emit: candidate_vcf
    tuple val(meta), path("*.candidateSmallIndels.vcf.gz"), path("*.candidateSmallIndels.vcf.gz.tbi") , emit: candidatesmallindels_vcf

    script:
    def prefix   = "${meta.id}.manta"
    def intervals_cmd  = intervals ? "--exome --callRegions ${intervals}" : ""
    """
    configManta.py \
        --bam ${bam} \
        --reference ${fasta} \
        ${intervals_cmd} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        ${prefix}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        ${prefix}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        ${prefix}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        ${prefix}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        ${prefix}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        ${prefix}.diploidSV.vcf.gz.tbi
    """
}


process MANTA_SOMATIC {
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
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path  fasta
    path  fai
    path  dict
    path intervals

    output:
    tuple val(meta), path("*.somaticSV.vcf.gz"), path("*.somaticSV.vcf.gz.tbi") , emit: somatic_vcf
    tuple val(meta), path("*.diploidSV.vcf.gz"), path("*.diploidSV.vcf.gz.tbi") , emit: diploid_vcf
    tuple val(meta), path("*.candidateSV.vcf.gz"), path("*.candidateSV.vcf.gz.tbi") , emit: candidate_vcf
    tuple val(meta), path("*.candidateSmallIndels.vcf.gz"), path("*.candidateSmallIndels.vcf.gz.tbi") , emit: candidatesmallindels_vcf

    script:
    def prefix   = "${meta.id}.manta"
    def intervals_cmd  = intervals ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    bgzip --threads ${task.cpus} -c ${intervals} > call_targets.bed.gz ; tabix call_targets.bed.gz

    configManta.py \
        --normalBam ${normal_bam} \
        --tumorBam ${tumor_bam} \
        --reference ${fasta} \
        ${intervals_cmd} \
        --runDir Manta
    python Manta/runWorkflow.py -m local -j ${task.cpus}
    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        ${prefix}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        ${prefix}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        ${prefix}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        ${prefix}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        ${prefix}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        ${prefix}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        ${prefix}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        ${prefix}.somaticSV.vcf.gz.tbi
    """
}
