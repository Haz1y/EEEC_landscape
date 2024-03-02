nextflow.enable.dsl = 2

include { create_fastq_channels; create_noidx_vcf_channels; create_merge_vcf_channels; create_somatic_bam_channels;create_varscan_vcf_channels; create_bam_channels;create_vcf_channels;create_absolute_channels; create_pyclone_channels } from '../modules/function'
include { VARSCAN } from '../modules/tools/varscan' addParams( outdir: params.outdir+"/1.varscan", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { VARSCAN_SOMATIC_FILTER } from '../modules/tools/varscan' addParams( outdir: params.outdir+"/2.varscan_filter", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { VARSCAN2MAF } from '../modules/tools/vcf2maf' addParams( outdir: params.outdir+"/3.varscan_filter_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { INDELOCATOR } from '../modules/tools/indelocator' addParams( outdir: params.outdir+"/3.indelocator", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { RESCUE_INDELOCATOR } from './python.nf' addParams( outdir: params.outdir+"/4.rescue_indelocator", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { INDELOCATORVCF2MAF } from '../modules/tools/vcf2maf' addParams( outdir: params.outdir+"/5.rescue_indelocator_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
//include { MANTA_SOMATIC } from '../modules/tools/manta' addParams( outdir: params.outdir+"/2.pindel", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { VCF_MERGE } from '../modules/tools/vcf_merge' addParams( outdir: params.outdir+"/11.merge_snv_indel_without_rescue_minN2_vcf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { COVERVCF2MAF } from '../modules/tools/vcf2maf' addParams( outdir: params.outdir+"/12.merge_snv_indel_without_rescue_minN2_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { UNFILTERINDELOCATORVCF2MAF } from '../modules/tools/vcf2maf' addParams( outdir: params.outdir+"/5.total_rescue_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)




workflow INDEL_CALLING {

  take:
  ch_varscan_vcf
 // ch_final_input
  //fasta
  //fasta_index
  //fasta_dict
  //intervals
  //gnomad
  //gnomad_idx
  //pon
  //pon_idx
  //variants_for_contamination
  //variants_for_contamination_idx
  //m2_extra_args

  main:

/*
  VARSCAN(
    ch_final_input,
    fasta,
    fasta_index,
    fasta_dict
    )

  INDELOCATOR(
    ch_final_input,
    fasta,
    fasta_index,
    fasta_dict
    )

    min_tumor_dp = 50
    min_normal_dp = 50
    min_tumor_af = 0.2
    max_normal_af = 0.05

  RESCUE_INDELOCATOR(
   ch_noidx_vcf,
   min_tumor_dp,
   min_normal_dp,
   min_tumor_af,
   max_normal_af
      )
*/

   VARSCAN_SOMATIC_FILTER(
    ch_varscan_vcf
        )

  emit:
  //  varscan_snv = VARSCAN.out.varscan_snv
  //  varscan_indel = VARSCAN.out.varscan_indel
  // indelocator_indel = INDELOCATOR.out.indelocator_indel
     varscan_somatic_hc  = VARSCAN_SOMATIC_FILTER.out.varscan_somatic_hc
}

workflow {

  genome = params.genome
  bwa_index = file(params.genomes[genome].bwa)
  fasta = file(params.genomes[genome].fasta)
  fasta_index = file(params.genomes[genome].fasta_fai)
  fasta_dict = file(params.genomes[genome].dict)
  intervals = params.exome ? file(params.genomes[genome].wes_target_bed) : file(params.genomes[genome].wgs_intervals)
  known_sites = file(params.genomes[genome].known_indels)
  known_sites_tbl = file(params.genomes[genome].known_indels_index)
  dbsnp = file(params.genomes[genome].dbsnp)
  dbsnp_index = file(params.genomes[genome].dbsnp_index)
  gnomad = file(params.genomes[genome].germline_resource)
  gnomad_idx = file(params.genomes[genome].germline_resource_index)
  variants_for_contamination = file(params.genomes[genome].variants_for_contamination)
  variants_for_contamination_idx = file(params.genomes[genome].variants_for_contamination_index)
  m2_extra_args = ""

  mode = params.mode

/*
  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_somatic_bam_channels(it) }
  .set { ch_bam }

  ch_bam.branch {
    meta, bam ->
    single  : bam.size() == 2
    return [ meta, bam[0], bam[1] ]
    paired  : bam.size() > 2
    return [ meta, bam[0], bam[1], bam[2], bam[3] ]
  }.set {ch_paired_bam}

  ch_paired_bam.paired
  .set{ch_final_input}



  INDEL_CALLING(
    ch_final_input,
    fasta,
    fasta_index,
    fasta_dict,
    intervals
    )
*/
    vcf2maf_path = "/data/person/huz/my_software/vcf2maf-1.6.21/"
    vep_path = "/data/person/huz/my_software/miniconda3/envs/DNAseq/bin"
    vep_data = "/data/person/huz/ref/vep/homo_sapiens/cache"

    ch_input = Channel.fromPath(params.input)
    ch_input
    .splitCsv( header:true, sep:',' )
    .map { create_noidx_vcf_channels(it) }
    .set { ch_noidx_vcf }

UNFILTERINDELOCATORVCF2MAF(
        ch_noidx_vcf,
        fasta,
        fasta_index,
        fasta_dict,
        vcf2maf_path,
        vep_path,
        vep_data
    )



/*
    min_tumor_dp = 50
    min_normal_dp = 50
    min_tumor_af = 0.2
    max_normal_af = 0.05


    INDEL_CALLING(
    ch_varscan_vcf
      )


    VARSCAN2MAF(
        INDEL_CALLING.out.varscan_somatic_hc,
        fasta,
        fasta_index,
        fasta_dict,
        vcf2maf_path,
        vep_path,
        vep_data
    )
    VCF_MERGE(
        ch_merge_vcf,
        fasta,
        fasta_index,
        fasta_dict
        )

    COVERVCF2MAF(
      VCF_MERGE.out.cover_vcf,
      fasta,
      fasta_index,
      fasta_dict,
      vcf2maf_path,
      vep_path,
      vep_data
      )
*/
}
