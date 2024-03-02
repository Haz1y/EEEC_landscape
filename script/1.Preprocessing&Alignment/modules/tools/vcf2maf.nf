
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process WY_VCF2MAF {
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
    tuple val(meta), path(vcf),path(tbi)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}"


    """
    bcftools view -f PASS ${vcf} > ${prefix}.mutect2_filtered.vcf
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${prefix}.mutect2_filtered.vcf \\
    --output-maf ${prefix}.unannovar.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${prefix}  --normal-id  ${prefix}_N --ncbi-build GRCh37

    rm -rf ${prefix}.mutect2_filtered.vcf
    """
}

process VCF2MAF {
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
    tuple val(meta), path(vcf),path(tbi)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}"
    def tumor      = "${meta.tumor}"
    def normal     = "${meta.normal}"

    """
    bcftools view -f PASS ${vcf} > ${prefix}.mutect2_filtered.vcf
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${prefix}.mutect2_filtered.vcf \\
    --output-maf ${prefix}.unannovar.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37

    rm -rf ${prefix}.mutect2_filtered.vcf
    """
}

process COVERVCF2MAF {
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
    tuple val(meta), path(cover_vcf)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}"
    def tumor      = "${meta.tumor}"
    def normal     = "${meta.normal}"

    """
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${cover_vcf} \\
    --output-maf ${prefix}.cover.minN2.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37
    """
}

process STRELKA2VCF2MAF_GRMC38 {
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
    tuple val(meta), path(snv_vcf), path(snv_tbi), path(indel_vcf),path(indel_tbi)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}"

    """
    bcftools view -f PASS ${snv_vcf} > ${prefix}.strelka2.pass.snv.vcf
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${prefix}.strelka2.pass.snv.vcf \\
    --output-maf ${prefix}.strelka2.pass.snv.maf  \\
    --ref-fasta ${fasta} \\
    --vcf-tumor-id TUMOR \\
    --vcf-normal-id NORMAL \\
    --species mus_musculus \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${prefix}  --normal-id  ${meta.patient} --ncbi-build GRCm38

    rm -rf ${prefix}.strelka2.pass.snv.vcf

    bcftools view -f PASS ${indel_vcf} > ${prefix}.strelka2.pass.indel.vcf
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${prefix}.strelka2.pass.indel.vcf \\
    --output-maf ${prefix}.strelka2.pass.indel.maf  \\
    --ref-fasta ${fasta} \\
    --vcf-tumor-id TUMOR \\
    --vcf-normal-id NORMAL \\
    --species mus_musculus \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${prefix}  --normal-id  ${meta.patient} --ncbi-build GRCm38

    rm -rf ${prefix}.strelka2.pass.indel.vcf
    """
}

process STRELKA2VCF2MAF {
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
    tuple val(meta), path(snv_vcf), path(snv_tbi), path(indel_vcf),path(indel_tbi)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}"

    """
    bcftools view -f PASS ${snv_vcf} > ${prefix}.strelka2.pass.snv.vcf
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${prefix}.strelka2.pass.snv.vcf \\
    --output-maf ${prefix}.strelka2.pass.snv.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${prefix}  --normal-id  ${meta.patient}_N --ncbi-build GRCh37

    rm -rf ${prefix}.strelka2.pass.snv.vcf

    bcftools view -f PASS ${indel_vcf} > ${prefix}.strelka2.pass.indel.vcf
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${prefix}.strelka2.pass.indel.vcf \\
    --output-maf ${prefix}.strelka2.pass.indel.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${prefix}  --normal-id  ${meta.patient}_N --ncbi-build GRCh37

    rm -rf ${prefix}.strelka2.pass.indel.vcf
    """
}

process VARSCAN2MAF {
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
    tuple val(meta), path(somatic_hc_snv_vcf), path(somatic_hc_indel_vcf)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}.varscan.somatic.hc"
    def tumor      = "${meta.tumor}"
    def normal     = "${meta.normal}"

    """
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${somatic_hc_snv_vcf} \\
    --output-maf ${prefix}.snv.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37

    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${somatic_hc_indel_vcf} \\
    --output-maf ${prefix}.indel.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37
    """
}

process VARSCAN2MAF_GRMC38 {
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
    tuple val(meta), path(somatic_hc_snv_vcf), path(somatic_hc_indel_vcf)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}.varscan.somatic.hc"


    """
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${somatic_hc_snv_vcf} \\
    --output-maf ${prefix}.snv.maf  \\
    --ref-fasta ${fasta} \\
    --vcf-tumor-id TUMOR \\
    --vcf-normal-id NORMAL \\
    --species mus_musculus \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${meta.id}  --normal-id  ${meta.patient} --ncbi-build GRCm38

    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${somatic_hc_indel_vcf} \\
    --output-maf ${prefix}.indel.maf  \\
    --ref-fasta ${fasta} \\
    --vcf-tumor-id TUMOR \\
    --vcf-normal-id NORMAL \\
    --species mus_musculus \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${meta.id}  --normal-id  ${meta.patient} --ncbi-build GRCm38
    """
}

process INDELOCATORVCF2MAF {
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
    tuple val(meta), path(somatic_only_vcf), path(rescue_only_vcf), path(passed_total_vcf)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}.indelocator"
    def tumor      = "${meta.tumor}"
    def normal     = "${meta.normal}"

    """
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${somatic_only_vcf} \\
    --output-maf ${prefix}.somatic_only.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37

    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${rescue_only_vcf} \\
    --output-maf ${prefix}.rescue_only.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37

    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${passed_total_vcf} \\
    --output-maf ${prefix}.passed_total.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37
    """
}

process UNFILTERINDELOCATORVCF2MAF {
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
    tuple val(meta), path(indelocator_vcf)
    path fasta
    path fasta_index
    path fasta_dict
    val vcf2maf_path
    val vep_path
    val vep_data

    output:
    tuple val(meta), path("*.maf"), emit: maf

    script:
    def prefix     = "${meta.id}.indelocator"
    def tumor      = "${meta.tumor}"
    def normal     = "${meta.normal}"

    """
    perl ${vcf2maf_path}/vcf2maf.pl \\
    --input-vcf ${indelocator_vcf} \\
    --output-maf ${prefix}.unfilter.maf  \\
    --ref-fasta ${fasta} \\
    --vep-path ${vep_path} \\
    --vep-data ${vep_data} \\
    --tumor-id ${tumor}  --normal-id  ${normal} --ncbi-build GRCh37
    """
}
