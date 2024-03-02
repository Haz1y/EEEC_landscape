
conda_env = (params.conda_env? params.conda_env : "bioconda::strelka=2.9.10")
singularity_url = "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
docker_url = "quay.io/biocontainers/strelka:2.9.10--0"

process VCF_MERGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/pyclone" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(mutect_vcf), path(strelka_snv), path(strelka_indel), path(varscan_snv), path(varscan_indel), path(rescue_indel)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("${meta.id}.cover.vcf"), emit: cover_vcf
    //更换format名字待定
    //varscan NORMAL  TUMOR
    //rescue TJ_L_026_C      TJ_L_026_N
    //mutect2 FD_AEH_001_C    FD_AEH_001_N
    //strelka2 NORMAL  TUMOR
    //less ${meta.id}.mutect2.pass.total.vcf | awk 'BEGIN{FS="\t";OFS="\t"}{if ($1 == "#CHROM") print $1, $2, $3, $4, $5, $6, $7, $8, $9, "TUMOR", "NORMAL"; else print}'  > ${meta.id}.mutect2.pass.total.rename.vcf
    //less ${rescue_indel}  | awk 'BEGIN{FS="\t";OFS="\t"}{if ($1 == "#CHROM") print $1, $2, $3, $4, $5, $6, $7, $8, $9, "TUMOR", "NORMAL"; else print}' > ${meta.id}.rescue.pass.total.rename.vcf

    script:
    def prefix   = "${meta.id}.cover"
    """
    /data/person/huz/my_software/miniconda3/envs/DNAseq/bin/bcftools view -f PASS ${mutect_vcf} > ${meta.id}.mutect2.pass.total.vcf
    /data/person/huz/my_software/miniconda3/envs/DNAseq/bin/bcftools view -f PASS ${strelka_snv} > ${meta.id}.strelka2.pass.snv.vcf
    /data/person/huz/my_software/miniconda3/envs/DNAseq/bin/bcftools view -f PASS ${strelka_indel} > ${meta.id}.strelka2.pass.indel.vcf

    sed -i s/${meta.id}/TUMOR/  ${meta.id}.mutect2.pass.total.vcf
    sed -i s/${meta.normal}/NORMAL/ ${meta.id}.mutect2.pass.total.vcf
    sed -i s/${meta.id}/TUMOR/ ${rescue_indel}
    sed -i s/${meta.normal}/NORMAL/ ${rescue_indel}

    GenomeAnalysisTK -R ${fasta} \\
    -T CombineVariants -o ${prefix}.vcf \\
    --variant:varscan ${varscan_snv} --variant:strelka ${meta.id}.strelka2.pass.snv.vcf \\
    --variant:varindel ${varscan_indel} --variant:sindel ${meta.id}.strelka2.pass.indel.vcf \\
    --variant:mutect ${meta.id}.mutect2.pass.total.vcf \\
    -minN 2 \\
    -genotypeMergeOptions PRIORITIZE -priority mutect,strelka,varscan,sindel,varindel
    """
}


/*
/data/person/huz/my_software/miniconda3/envs/DNAseq/bin/bcftools view -f PASS ${mutect_vcf} > ${meta.id}.mutect2.pass.total.vcf
/data/person/huz/my_software/miniconda3/envs/DNAseq/bin/bcftools view -f PASS ${strelka_snv} > ${meta.id}.strelka2.pass.snv.vcf
/data/person/huz/my_software/miniconda3/envs/DNAseq/bin/bcftools view -f PASS ${strelka_indel} > ${meta.id}.strelka2.pass.indel.vcf

sed -i s/${meta.id}/TUMOR/  ${meta.id}.mutect2.pass.total.vcf
sed -i s/${meta.normal}/NORMAL/ ${meta.id}.mutect2.pass.total.vcf
sed -i s/${meta.id}/TUMOR/ ${rescue_indel}
sed -i s/${meta.normal}/NORMAL/ ${rescue_indel}

GenomeAnalysisTK -R ${fasta} \\
-T CombineVariants -o ${prefix}.vcf \\
--variant:varscan ${varscan_snv} --variant:strelka ${meta.id}.strelka2.pass.snv.vcf \\
--variant:varindel ${varscan_indel} --variant:sindel ${meta.id}.strelka2.pass.indel.vcf \\
--variant:mutect ${meta.id}.mutect2.pass.total.vcf --variant:rescue ${rescue_indel} \\
-minN 2 \\
-genotypeMergeOptions PRIORITIZE -priority mutect,strelka,varscan,sindel,rescue,varindel
*/
