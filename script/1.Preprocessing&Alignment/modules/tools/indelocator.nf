
conda_env = (params.conda_env? params.conda_env : "bioconda::strelka=2.9.10")
singularity_url = "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
docker_url = "quay.io/biocontainers/strelka:2.9.10--0"

process INDELOCATOR {
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
    tuple val(meta), path("${meta.id}.indelocator.indels.vcf"), emit: indelocator_indel

    script:
    def prefix   = "${meta.id}.indelocator"
    """
    /usr/lib/jvm/java-1.7.0-openjdk-1.7.0.171-2.6.13.2.el7.x86_64/jre/bin/java -jar \
      /data/person/huz/my_software/Indelocater/IndelGenotyper.36.3336-GenomeAnalysisTK.jar \
      -T IndelGenotyperV2 -I:normal ${normal_bam} -I:tumor ${tumor_bam} \
      -R ${fasta} \\
      -somatic -minCoverage 6 -minNormalCoverage 4 -minFraction 0.2 \\
      -o ${prefix}.indels.vcf
    """
}
