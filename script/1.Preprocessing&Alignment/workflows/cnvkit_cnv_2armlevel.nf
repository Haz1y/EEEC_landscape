nextflow.enable.dsl = 2

include { create_fastq_channels; create_noidx_vcf_channels; create_seg_channels; create_merge_vcf_channels; create_somatic_bam_channels;create_varscan_vcf_channels; create_bam_channels;create_vcf_channels;create_absolute_channels; create_pyclone_channels } from '../modules/function'
include { CNVKIT_ASCETS } from '../modules/tools/ascets' addParams( outdir: params.outdir+"/25_cnvkit_ascets", publish_enable: true, publish_dir_mode: params.publish_dir_mode)





workflow TO_DO {

  take:
  ch_seg
  genomic_arm_coordinates
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

  ASCETS(
      ch_seg,
      genomic_arm_coordinates
      )
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
  VARSCAN_SOMATIC_FILTER(
   ch_varscan_vcf
       )
*/



  //emit:
  //  varscan_snv = VARSCAN.out.varscan_snv
  //  varscan_indel = VARSCAN.out.varscan_indel
  // indelocator_indel = INDELOCATOR.out.indelocator_indel
    // varscan_somatic_hc  = VARSCAN_SOMATIC_FILTER.out.varscan_somatic_hc
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


    ch_input = Channel.fromPath(params.input)
    ch_input
    .splitCsv( header:true, sep:',' )
    .map { create_seg_channels(it) }
    .set { ch_seg }



    genomic_arm_coordinates = "/data/person/huz/my_software/ascets/genomic_arm_coordinates_hg19.txt"
    remove_bash = "/data/person/huz/DNAseq/CNVkit_UCEC/06_somatic_cnv_calling/removeY.bash"

    CNVKIT_ASCETS(
    ch_seg,
    genomic_arm_coordinates

    )

}
