nextflow.enable.dsl = 2

include { create_somatic_bam_channels } from '../modules/function'

include { GATK4_MUTECT2} from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { GATK4_CALCULATECONTAMINATION } from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_LEARNREADORIENTATIONMODEL } from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { GATK4_FILTERMUTECTCALLS } from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)


workflow MUTECT2_PAIRED {
    /*
    ch_bam  meta [id, sample, type, ], tumor_bam, tumor_bai, normal_bam, normal_bai
    */
  take:
  ch_bam
  fasta
  fasta_index
  fasta_dict
  intervals
  gnomad
  gnomad_idx
  variants_for_contamination
  variants_for_contamination_idx
  m2_extra_args

  main:

  GATK4_MUTECT2(
      ch_bam,
      fasta,
      fasta_index,
      fasta_dict,
      intervals,
      gnomad,
      gnomad_idx,
      variants_for_contamination,
      variants_for_contamination_idx,
      m2_extra_args
      )

  // TODO TUMOR_ONLY
  //GATK4_MUTECT2.out.vcf.concat(GATK4_MUTECT2_TUMOR_ONLY.out.vcf)
  //.set{ ch_mutect_vcf }

  //GATK4_MUTECT2.out.stats.concat(GATK4_MUTECT2_TUMOR_ONLY.out.stats)
  //.set{ ch_mutect_stats }

  //GATK4_MUTECT2_TUMOR_ONLY.out.pileup.map {
  //  meta, pileup ->
  //  [meta, pileup, null]
  //}
  //.mix(GATK4_MUTECT2.out.pileup)
  //.set{ ch_mutect_pileup }

  ch_mutect_pileup = GATK4_MUTECT2.out.pileup
  ch_mutect_stats  = GATK4_MUTECT2.out.stats
  ch_mutect_vcf    = GATK4_MUTECT2.out.vcf
  ch_mutect_f1r2   = GATK4_MUTECT2.out.f1r2

  GATK4_CALCULATECONTAMINATION(
    ch_mutect_pileup
    )

  GATK4_LEARNREADORIENTATIONMODEL(
    ch_mutect_f1r2
    )

  ch_mutect_vcf.concat(
    GATK4_CALCULATECONTAMINATION.out.contamination,
    GATK4_CALCULATECONTAMINATION.out.maf_segments,
    ch_mutect_stats,
    GATK4_LEARNREADORIENTATIONMODEL.out.artifact_prior_table
    )
    .groupTuple(by: [0])
    .map {
      meta, v1, v2 ->
      [meta, v1[0], v2, v1[1], v1[2], v1[3], v1[4]]
      }
    .set { ch_filter }

  GATK4_FILTERMUTECTCALLS(
    ch_filter,
    fasta,
    fasta_index,
    fasta_dict
    )

  emit:
    vcf = GATK4_FILTERMUTECTCALLS.out.vcf
  //  DP >=8 VAF >= 0.05 AD >=5 tumor
  //  VAF < 0.01  NAT   strand bias <= 0.95


}


workflow {

  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_somatic_bam_channels(it) }
  .set { ch_bam }

  genome = params.genome
  fasta = file(params.genomes[genome].fasta)
  fasta_index = file(params.genomes[genome].fasta_fai)
  fasta_dict = file(params.genomes[genome].dict)
  intervals = params.exome ? file(params.genomes[genome].wes_target_bed) : file(params.genomes[genome].wgs_intervals)
  gnomad = file(params.genomes[genome].germline_resource)
  gnomad_idx = file(params.genomes[genome].germline_resource_index)
  variants_for_contamination = file(params.genomes[genome].variants_for_contamination)
  variants_for_contamination_idx = file(params.genomes[genome].variants_for_contamination_index)
  m2_extra_args = ""

  ch_bam
  .branch {
    meta, bam ->
    single  : bam.size() == 2
    return [ meta, bam[0], bam[1]]
    paired  : bam.size() > 2
    return [ meta, bam[0], bam[1], bam[2], bam[3]]
  }
  .set { ch_paired_bam }

  MUTECT2_PAIRED(
    ch_paired_bam.paired,
    fasta,
    fasta_index,
    fasta_dict,
    intervals,
    gnomad,
    gnomad_idx,
    variants_for_contamination,
    variants_for_contamination_idx,
    m2_extra_args
    )
}
