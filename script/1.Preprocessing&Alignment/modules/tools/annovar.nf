
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process ANNOVAR {
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
    tuple val(meta), path(vcf), path(vcf_index)
    path fasta
    path fasta_index
    path fasta_dict
    val genome
    val annovar_path

    output:
    tuple val(meta), path("*.annovar.hg19anno.csv"), emit: csv

    script:
    def prefix     = "${meta.id}.annovar"

    """

    bcftools norm -m-both -o ex1.step1.vcf ${vcf}
    bcftools norm -f ${fasta} -o ex1.step2.vcf ex1.step1.vcf
    ${annovar_path}/convert2annovar.pl -format vcf4 -filter pass ex1.step2.vcf > ex2.avinput
    ${annovar_path}/table_annovar.pl -protocol refGene --operation g -nastring . -csvout -buildver  ${genome} ex2.avinput ${annovar_path}/humandb/
    mv ex2.avinput.${genome}_multianno.csv ${prefix}.${genome}anno.csv

    """
}

process GERMLINE_ANNOVAR {
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
    tuple val(meta), path(snp_vcf), path(indel_vcf), path(raw_vcf), path(raw_tbi)
    path fasta
    path fasta_index
    path fasta_dict
    val annovar_path

    output:
    tuple val(meta), path("*.hg19_multianno.txt"), emit: annovar_germ_tsv

    script:
    def prefix     = "${meta.id}.annovar"

    """
    ${annovar_path}/convert2annovar.pl -format vcf4 \\
      ${snp_vcf} \\
      -outfile ${meta.id}.snps.annovar_raw.tsv \\
      -include \\
      -withzyg

    ${annovar_path}/table_annovar.pl ${meta.id}.snps.annovar_raw.tsv \\
      /data/person/huz/ref/annovar/humandb \\
      -buildver hg19 \\
      -out ${meta.id}.snps.vqsr.annovar.tsv\\
      -remove \\
      -protocol refGene,avsnp150,clinvar_20220320 \\
      -operation g,f,f \\
      -nastring NAN

    ${annovar_path}/convert2annovar.pl -format vcf4 \\
      ${indel_vcf} \\
      -outfile ${meta.id}.indels \\
      -include \\
      -withzyg

    ${annovar_path}/table_annovar.pl ${meta.id}.indels.annovar_raw.tsv \\
      /data/person/huz/ref/annovar/humandb \\
      -buildver hg19 \\
      -out ${meta.id}.indels.vqsr \\
      -remove \\
      -protocol refGene,avsnp150,clinvar_20220320 \\
      -operation g,f,f \\
      -nastring NAN

      ${annovar_path}/convert2annovar.pl -format vcf4 \\
        ${raw_vcf} \\
        -outfile ${meta.id}.raw_vcf.annovar_raw.tsv \\
        -include \\
        -withzyg

      ${annovar_path}/table_annovar.pl ${meta.id}.raw_vcf.annovar_raw.tsv \\
        /data/person/huz/ref/annovar/humandb \\
        -buildver hg19 \\
        -out ${meta.id}.raw_vcf \\
        -remove \\
        -protocol refGene,avsnp150,clinvar_20220320 \\
        -operation g,f,f \\
        -nastring NAN

    """
}

process GERMLINE_DIRECT_ANNOVAR {
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
    tuple val(meta), path(raw_vcf), path(raw_tbi)
    path fasta
    path fasta_index
    path fasta_dict
    val annovar_path

    output:
    tuple val(meta), path("*.hg19_multianno.txt"), emit: annovar_germ_tsv

    script:
    def prefix     = "${meta.id}.annovar"

    """
      ${annovar_path}/convert2annovar.pl -format vcf4 \\
        ${raw_vcf} \\
        -outfile ${meta.id}.raw_vcf.annovar_raw.tsv \\
        -include \\
        -withzyg

      ${annovar_path}/table_annovar.pl ${meta.id}.raw_vcf.annovar_raw.tsv \\
        /data/person/huz/ref/annovar/humandb \\
        -buildver hg19 \\
        -out ${meta.id}.raw_vcf \\
        -remove \\
        -protocol refGene,avsnp150,clinvar_20220320 \\
        -operation g,f,f \\
        -nastring NAN

    """
}
