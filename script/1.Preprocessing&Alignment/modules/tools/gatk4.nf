
conda_env = (params.conda_env? params.conda_env: "bioconda::gatk4=4.2.0.0")
singularity_url = "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
docker_url = "quay.io/biocontainers/gatk4:4.2.0.0--0"

process GATK4_MARKDUPLICATES {
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
    tuple val(meta), path(bam)
    val   extra_args

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bai")    , emit: bam
    tuple val(meta), path("*.markdup.metrics"), emit: metrics

    script:
    def prefix   = "${meta.id}"
    def extra_args_command = extra_args ? "${extra_args}" : ""
    """
    gatk MarkDuplicates \\
        --INPUT ${bam} \\
        --METRICS_FILE ${prefix}.markdup.metrics \\
        --TMP_DIR . \\
        --ASSUME_SORT_ORDER coordinate \\
        --CREATE_INDEX true \\
        --OUTPUT ${prefix}.markdup.bam \\
        $extra_args

    """
}


process GATK4_APPLYBQSR {
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
    tuple val(meta), path(bam), path(bai), path(bqsr_table)
    path  fasta
    path  fastaidx
    path  dict
    path  intervals

    output:
    tuple val(meta), path("*.bqsr.bam"), path("*.bqsr.bai"), emit: bam

    script:
    def prefix   = "${meta.id}"
    def interval = intervals ? "-L ${intervals}" : ""
    """
    gatk ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${bqsr_table} \\
        --create-output-bam-index true \\
        ${interval} \\
        -O ${prefix}.bqsr.bam
    """
}


process GATK4_BASERECALIBRATOR {
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
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path intervalsBed
    path knownsites
    path knownsites_tbi
    path dbsnp
    path dbsnp_index

    output:
    tuple val(meta), path("*.recal.table"), emit: table

    script:
    def prefix   = "${meta.id}"
    def intervals_command = intervalsBed ? "-L ${intervalsBed}" : ""
    def sites_command = knownsites.collect{"--known-sites ${it}"}.join(' ')
    """
    gatk BaseRecalibrator  \\
        -R ${fasta} \\
        -I ${bam} \\
        --known-sites ${dbsnp} \\
        ${sites_command} \\
        ${intervals_command} \\
        --tmp-dir . \\
        -O ${prefix}.recal.table

    """


}


process GATK4_CALCULATECONTAMINATION {
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
    tuple val(meta), path(tumor_pileup), file(normal_pileup)

    output:
    tuple val(meta), path("*.contamination.table"), emit: contamination
    tuple val(meta), path("*.tumor_segment.table")   , emit: maf_segments

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CalculateContamination \\
        -I ${tumor_pileup} \\
        -matched ${normal_pileup} \\
        -O ${prefix}.contamination.table \\
        --tumor-segmentation ${prefix}.tumor_segment.table

    """
}

/*  TODO
process GATK4_FILTERALIGNMENTARTIFACTS {
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
  tuple val(meta), path(tumor_bam), path(tumor_bai), path(realignment_index_bundle), path(vcf), path(vcf_tbi)
  path fasta
  path fai
  path dict

  output:
  tuple val(meta), path("*.mutect2_filtered.vcf.gz"), emit: vcf
  tuple val(meta), path("*.mutect2_filtered.stats"), emit: stats

  script:
  def prefix    = "${meta.id}.mutect2_filtered"

  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
      FilterAlignmentArtifacts \\
      -R ${fasta} \\
      -V ${vcf} \\
      -I ${tumor_bam} \\
      --bwa-mem-index-image ${realignment_index_bundle} \
            ~{realignment_extra_args} \
            -O ~{output_vcf} \
      -O ${prefix}.vcf.gz
  """

}
*/

process GATK4_FILTERMUTECTCALLS {
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
    tuple val(meta), path(unfiltered_vcf), path(unfiltered_vcf_idx), path(contamination_table), path(maf_segments), path(mutect_stats), path(artifact_priors)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.mutect2_filtered.vcf.gz"), path("*.mutect2_filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.mutect2_filtered.stats"), emit: stats

    script:
    def prefix    = "${meta.id}.mutect2_filtered"

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        FilterMutectCalls \\
        -R ${fasta} \\
        -V ${unfiltered_vcf} \\
        --contamination-table ${contamination_table} \\
        --tumor-segmentation ${maf_segments} \\
        --stats ${mutect_stats} \\
        --ob-priors ${artifact_priors} \\
        --filtering-stats  ${prefix}.stats \\
        -O ${prefix}.vcf.gz
    """
}


process GATK4_HAPLOTYPECALLER {
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
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path intervals
    path dbsnp
    path dbsnp_tbl
    val use_gvcf

    output:
    tuple val(meta), path("*.vcf.gz"),path("*.tbi") , emit: vcf

    script:
    def prefix    = "${meta.id}"
    def mode = use_gvcf ? "-ERC GVCF" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -L ${intervals} \\
        ${mode} \\
        -O ${prefix}.vcf.gz \\
        -D ${dbsnp}
    """
}

process GATK4_GENOMICSDBIIMPORT {
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
    tuple val(meta), path(gvcfs), path(tbis)
    path fasta
    path fai
    path dict
    path intervals

    output:
    tuple val(meta), path("gvcfs_db/"), emit: vcf_db

    script:
    def prefix    = "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        GenomicsDBImport \\
        -R ${fasta} \\
        ${gvcfs.flatten().collect{ "-V $it " }.join()} \\
        -L ${intervals} \\
        --merge-input-intervals TRUE \\
        --genomicsdb-workspace-path gvcfs_db
    """
}

process GATK4_GENOTYPEGVCFS {
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
    tuple val(meta), path(vcf_db)
    path fasta
    path fai
    path dict
    path intervals
    path dbsnp
    path dbsnp_index

    output:
    tuple val(meta), path("germline_merge.vcf.gz"), emit: merge_vcf

    script:
    def prefix    = "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        GenotypeGVCFs \\
        -R ${fasta} \\
        -V gendb://gvcfs_db \\
        -D ${dbsnp} \\
        -L ${intervals} \\
        -O germline_merge.vcf.gz
    """
}

process GATK4_TOTAL_GERMLINE_VQSR {
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
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fai
    path dict
    path known_indels
    path known_indels_index

    output:
    tuple val(meta), path("${meta.id}.snps.VQSR.filter.vcf"), path("${meta.id}.indels.VQSR.filter.vcf"),path("${vcf}"), path("${tbi}"),emit: vqsr_vcf
    tuple val(meta), path("${meta.id}.snps.plots.R"), path("${meta.id}.indels.plots.R"),emit: vqsr_plot_r

    script:
    def prefix   = "${meta.id}"
    """
    gatk VariantRecalibrator  \\
        -R ${fasta} \\
        -V ${vcf} \\
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/person/huz/ref/gatk/hapmap_3.3.b37.vcf.gz \\
        -resource:omini,known=false,training=true,truth=false,prior=12.0 /data/person/huz/ref/gatk/1000G_omni2.5.b37.vcf.gz \\
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/person/huz/ref/gatk/1000G_phase1.snps.high_confidence.b37.vcf.gz \\
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/person/huz/ref/gatk/dbsnp_138.b37.vcf.gz \\
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
        -mode SNP \\
        --max-gaussians 3 \\
        --rscript-file ${prefix}.snps.plots.R \\
        --tranches-file ${prefix}.snps.tranches \\
        -O ${prefix}.snps.recal && \\
   gatk ApplyVQSR  \\
       -R ${fasta} \\
       -V ${vcf} \\
       --ts_filter_level 99.0 \\
       --tranches-file ${prefix}.snps.tranches \\
       --recal-file ${prefix}.snps.recal \\
       -mode SNP \\
       -O ${prefix}.snps.VQSR.vcf.gz

   gatk VariantRecalibrator  \\
       -R ${fasta} \\
       -V ${vcf} \\
       -resource:mills,known=true,training=true,truth=true,prior=12.0 /data/database/homo_sapiens/GRCh37/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf \\
       -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/database/homo_sapiens/GRCh37/Annotation/GATKBundle/dbsnp_138.b37.vcf \\
       -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
       -mode INDEL \\
       --max-gaussians 3 \\
       --rscript-file ${prefix}.indels.plots.R \\
       --tranches-file ${prefix}.indels.tranches \\
       -O ${prefix}.indels.recal && \\
  gatk ApplyVQSR  \\
      -R ${fasta} \\
      -V ${vcf} \\
      -ts-filter-level 99.0 \\
      --tranches-file ${prefix}.indels.tranches \\
      --recal-file ${prefix}.indels.recal \\
      -mode INDEL \\
      -O ${prefix}.indels.VQSR.vcf.gz

    bcftools view -f PASS ${prefix}.snps.VQSR.vcf.gz > ${prefix}.snps.VQSR.filter.vcf
    bcftools view -f PASS ${prefix}.indels.VQSR.vcf.gz > ${prefix}.indels.VQSR.filter.vcf

    """


}

process GATK4_GERMLINE_HARD_FILTER {
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
    tuple val(meta), path(vcf), path(tbi)
    val mode
    path fasta
    path fai
    path dict
    path intervalsBed
    path knownsites
    path knownsites_tbi
    path dbsnp
    path dbsnp_index

    output:
    tuple val(meta), path("*.recal.table"), emit: table

    script:
    def prefix   = "${meta.id}"
    def intervals_command = intervalsBed ? "-L ${intervalsBed}" : ""
    def sites_command = knownsites.collect{"--known-sites ${it}"}.join(' ')
    """
    gatk VariantFiltration  \\
        -R ${fasta} \\
        -V ${snp_vcf} \\
        --output ${prefix}_snp.filter.vcf \\
        --filter-name "SNP_DQ" \\
        --filter-expression "DQ < 2.0" \\
        --filter-name "SNP_MQ" \\
        --filter-expression "MQ < 40.0" \\
        --filter-name "SNP_FS" \\
        --filter-expression "FS > 60.0" \\
        --filter-name "SNP_SOR" \\
        --filter-expression "SOR > 3.0" \\
        --filter-name "SNP_MQRankSum" \\
        --filter-expression "MQRankSum < -12.5" \\
        --filter-name "SNP_ReadPosRankSum" \\
        --filter-expression "ReadPosRankSum < -8.0"

    gatk VariantFiltration  \\
        -R ${fasta} \\
        -V ${indel_vcf} \\
        --output ${prefix}_indel.filter.vcf \\
        --filter-name "INDEL_DQ" \\
        --filter-expression "DQ < 2.0" \\
        --filter-name "INDEL_FS" \\
        --filter-expression "FS > 200.0" \\
        --filter-name "INDEL_SOR \\
        --filter-expression "SOR > 10.0" \\
        --filter-name "INDEL_ReadPosRankSum" \\
        --filter-expression "ReadPosRankSum < -20.0" \\
        --filter-name "INDEL_InbreedingCoeff" \\
        --filter-expression "InbreedingCoeff < -0.8"
    """


}

process GATK4_VARIANTCALIBRATOR {
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
    tuple val(meta), path(vcf), path(tbi)
    val mode
    path fasta
    path fai
    path dict
    path intervalsBed
    path knownsites
    path knownsites_tbi
    path dbsnp
    path dbsnp_index

    output:
    tuple val(meta), path("*.recal.table"), emit: table

    script:
    def prefix   = "${meta.id}"
    def intervals_command = intervalsBed ? "-L ${intervalsBed}" : ""
    def sites_command = knownsites.collect{"--known-sites ${it}"}.join(' ')
    """
    gatk VariantRecalibrator  \\
        -R ${fasta} \\
        -V ${vcf} \\
        --output ${prefix}.indel.${mode}.reccal \\
        --tranches-file ${prefix}.${mode}.tranches \\
        --rscript-file ${prefix}.${mode}.plots.R \\
        --trust-all-polymorphic \\
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode ${mode} \\
        ${sites_command}
    """


}

process GATK4_LEARNREADORIENTATIONMODEL {

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
  tuple val(meta), path(f1r2)

  output:
  tuple val(meta), path("*.artifact-priors.tar.gz"), emit: artifact_prior_table

  script:
  def prefix    = "${meta.id}"

  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
    LearnReadOrientationModel \\
    -I ${f1r2} \\
    -O ${prefix}.artifact-priors.tar.gz

  """
}


process GATK4_MUTECT2 {
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
    path fai
    path dict
    path intervals
    path gnomad
    path gnomad_index
    path variants_for_contamination
    path variants_for_contamination_idx
    val extra_args

    output:
    tuple val(meta), path("*.vcf"), path("*.vcf.idx"), emit: vcf
    tuple val(meta), path("*.vcf.stats"), emit: stats
    tuple val(meta), path("*.tumor_pileup.table"), path("*.normal_pileup.table"), emit: pileup
    tuple val(meta), path("*.f1r2.tar.gz"), emit: f1r2

    script:

    def prefix    = "${meta.id}.mutect2"
    def pon_command = params.pon ? "--panel-of-normals ${pon}" : ""
    //def intervals_command = intervals ? "-L ${intervals}" : ""

    """
    gatk GetSampleName -R ${fasta} -I ${tumor_bam} -O tumor_name.txt
    tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"
    gatk GetSampleName -R ${fasta} -I ${normal_bam} -O normal_bam.txt
    normal_command_line="-I ${normal_bam} -normal `cat normal_bam.txt`"

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
      Mutect2 \\
      -R ${fasta} \\
      \${tumor_command_line} \\
      \${normal_command_line} \\
      -L ${intervals} \\
      --germline-resource ${gnomad} \\
      ${pon_command} \\
      --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
      ${extra_args} \\
      -O ${prefix}.vcf

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
      GetPileupSummaries \\
      -R ${fasta}\\
      -I ${tumor_bam} \\
      -L ${variants_for_contamination}  \\
      -V ${variants_for_contamination} \\
      -O ${prefix}.tumor_pileup.table

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        GetPileupSummaries \\
        -R ${fasta}\\
        -I ${normal_bam} \\
        -L ${variants_for_contamination}  \\
        -V ${variants_for_contamination} \\
        -O ${prefix}.normal_pileup.table
    """

}


process GATK4_MUTECT2_TUMOR_ONLY {
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
    tuple val(meta), path(tumor_bam), path(tumor_bai)
    path fasta
    path fai
    path dict
    path intervals
    path gnomad
    path gnomad_tbi
    file(pon)
    file(pon_tbi)

    output:
    tuple val(meta), path("*.vcf"), path("*.vcf.idx"), emit: vcf
    tuple val(meta), path("*.vcf.stats"), emit: stats
    tuple val(meta), path("*.tumor_pileups.table"), emit: pileup

    script:
    def prefix    = "${meta.id}.mutect2"
    def pon_command = params.pon ? "--panel-of-normals ${pon}" : ""
    def intervals_command = intervals ? "-L ${intervals}" : ""

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
      Mutect2 \\
      -R ${fasta}\\
      -I ${tumor_bam} \\
      ${intervals_command} \\
      --germline-resource ${gnomad} \\
      ${pon_command} \\
      -O ${prefix}.vcf

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
      GetPileupSummaries \\
      -R ${fasta}\\
      -I ${tumor_bam} \\
      ${intervals_command} \\
      -V ${gnomad} \\
      -O ${prefix}.tumor_pileups.table

    """
}


process GATK4_SELECTVARIANTS  {
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
    tuple val(meta), path(vcf), path(vcf_tbi)
    val extra_args

    output:
    tuple val(meta), path("*.vcf"), path("*.vcf.idx"), emit: vcf

    script:
    def prefix    = "${meta.id}"

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
      SelectVariants  \\
      -V ${vcf} \\
      ${extra_args} \\
      -O ${prefix}.vcf

    #gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
    #  IndexFeatureFile \\
    #  -R ${fasta}\\
    #  -I ${tumor_bam} \\
    #  ${intervals_command} \\
    #  -V ${gnomad} \\
    #  -O ${prefix}.tumor_pileups.table

    """
}
