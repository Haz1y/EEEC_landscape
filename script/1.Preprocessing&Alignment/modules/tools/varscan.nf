
conda_env = (params.conda_env? params.conda_env : "bioconda::strelka=2.9.10")
singularity_url = "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
docker_url = "quay.io/biocontainers/strelka:2.9.10--0"

process VARSCAN {
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

    output:
    tuple val(meta), path("${meta.id}.varscan.snv.vcf"),path("${meta.id}.varscan.indel.vcf"), emit: varscan_vcf

    script:
    def prefix   = "${meta.id}.pindel"
    def exome_cmd  = params.exome ? "--exome" : ""
    """
    samtools mpileup \\
    -q 1 \\
    -f ${fasta} \\
    ${normal_bam} ${tumor_bam} > ${meta.id}.normal_tumor.mpileup

    varscan somatic -Xmx${task.memory.toGiga()}g \\
    ${meta.id}.normal_tumor.mpileup \\
    --mpileup 1 \\
    --min-coverage-normal 8 --min-coverage-tumor 14 \\
    --min-var-freq 0.05  \\
    --strand-filter 1 \\
    --output-vcf 1 \\
    --output-snp ${meta.id}.varscan.snv.vcf \\
    --output-indel ${meta.id}.varscan.indel.vcf

    rm -rf ${meta.id}.normal_tumor.mpileup
    """
}

process VARSCAN_SOMATIC_FILTER {
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
    tuple val(meta), path(snv_vcf), path(indel_vcf)


    output:
    tuple val(meta), path("*.snv.Somatic.hc.vcf"), path("*.indel.Somatic.hc.vcf"), emit: varscan_somatic_hc
    tuple val(meta), path("*.snv.Somatic.vcf"), path("*.indel.Somatic.vcf"), emit: varscan_somatic_total
    tuple val(meta), path("*.Germline.*"), emit: varscan_germline_total
    tuple val(meta), path("*.LOH.*"),  emit: varscan_loh_total

    script:

    """
    varscan processSomatic ${snv_vcf} \\
      --min-tumor-freq  0.05   --max-normal-freq  0.05   --p-value  0.05

    varscan processSomatic ${indel_vcf} \\
    --min-tumor-freq  0.05   --max-normal-freq  0.05   --p-value  0.05
    """
}
