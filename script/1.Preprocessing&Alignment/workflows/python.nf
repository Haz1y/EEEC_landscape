conda_env = (params.conda_env? params.conda_env: "pysam python=3.6")
singularity_url = ""
docker_url = ""

process FILTER_MUTEC2_VCF {
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
    val   min_tumor_dp
    val   min_normal_dp
    val   min_tumor_ad
    val   min_tumor_af
    val   max_normal_af

    output:
    tuple val(meta), path("*.mutect2_filtered.vcf"), path("*.mutect2_filtered.vcf.idx"), emit: vcf

    script:
    def prefix     = "${meta.id}.annovar"

    """
    filter_paired_mutect_vcf.py --filter PASS --min_tumor_dp ${min_tumor_dp} \\
        --min_normal_dp ${min_normal_dp} --min_tumor_ad ${min_tumor_ad} --min_tumor_af ${min_tumor_af} --max_normal_af ${max_normal_af} \\
        -o ${prefix}.mutect2_filtered.vcf ${vcf}
    gatk IndexFeatureFile -I ${prefix}.mutect2_filtered.vcf
    """
}

process CREATE_SOMATIC_BAM_TABLE {
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
    path input_csv

    output:
    path("somatic_bam.csv"), emit: table

    script:

    """
    create_somatic_bam_table.py ${input_csv} somatic_bam.csv
    """
}

process RESCUE_INDELOCATOR {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta), path(vcf)
    val   min_tumor_dp
    val   min_normal_dp
    val   min_tumor_af
    val   max_normal_af

    output:
    tuple val(meta), path("*.somatic_only.vcf"), path("*.rescue_only.vcf"), path("*.total.vcf"), emit: rescue


    script:

    """
    python /data/person/huz/pipeline/nf-wgs/bin/rescue_indelocator.py ${vcf} --min_tumor_dp ${min_tumor_dp} \\
        --min_normal_dp ${min_normal_dp} --min_tumor_af ${min_tumor_af} --max_normal_af ${max_normal_af} \\
        --tumor_ID ${meta.id}
    """
}
