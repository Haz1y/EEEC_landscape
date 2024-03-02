
conda_env = (params.conda_env? params.conda_env : "bioconda::strelka=2.9.10")
singularity_url = "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
docker_url = "quay.io/biocontainers/strelka:2.9.10--0"

process PINDEL {
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
    tuple val(meta), path("${prefix}.vcf"), emit: pindel_vcf

    script:
    def prefix   = "${meta.id}.pindel"
    def exome_cmd  = params.exome ? "--exome" : ""
    """
    echo -e "${tumor_bam}\t150\t${meta.id}" > pindel.config
    echo -e "${normal_bam}\t150\t${meta.patient}_N" >> pindel.config

    pindel -T ${task.cpus} -f ${fasta} -i pindel.config \\
    -o ${meta.id}.pindel_raw -m 6

    pindel2vcf -p  ${meta.id}.pindel_raw -r ${fasta} -R 1000GenomesPilot-NCBI37 \\
             -d 20220323 -v ${prefix}.vcf -G

    """
}
