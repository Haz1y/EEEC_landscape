nextflow.enable.dsl = 2

include { create_somatic_bam_channels;create_OV_somatic_bam_channels } from '../modules/function'

include { GATK4_COLLECTREADCOUNTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir, publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { GATK4_CREATEREADCOUNTPANELOFNORMALS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_COLLECTALLELICCOUNTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir+"/01_alleic_counts", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_DENOISEREADCOUNTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir+"/02_denoised_counts", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_MODELSEGMENTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_CALLCOPYRATIOSEGMENTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir+"/04_called_seg", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_POLTDENOISEDRATIOS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_POLTMODELSEGMENTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_FUNCOCATESEGMENTS} from '../modules/tools/gatk_cnv' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)


workflow SOMATIC_CNV {
    /*
    ch_bam  meta [id, sample, type, ], tumor_bam, tumor_bai, normal_bam, normal_bai
    */
  take:
  ch_bam
  fasta
  fasta_index
  fasta_dict
  process_intervals
  common_sites

  main:

  GATK4_COLLECTREADCOUNTS(
      ch_bam,
      fasta,
      fasta_index,
      fasta_dict,
      process_intervals
      )
  ch_tumor_counts = GATK4_COLLECTREADCOUNTS.out.tumor_hdf5
  ch_normal_counts = GATK4_COLLECTREADCOUNTS.out.normal_hdf5

  GATK4_COLLECTREADCOUNTS.out.tumor_hdf5.concat(
    GATK4_COLLECTREADCOUNTS.out.normal_hdf5
    )
    .groupTuple(by: [0])
    .map {
      meta, v1 ->
      [meta, v1[0], v1[1] ]
      }
    .set { ch_counts }

  GATK4_CREATEREADCOUNTPANELOFNORMALS(
      ch_normal_counts.collect{it[1]}
      )
  pon = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.pon

  GATK4_COLLECTALLELICCOUNTS(
      ch_bam,
      fasta,
      fasta_index,
      fasta_dict,
      common_sites
      )
  ch_tumor_allelic = GATK4_COLLECTALLELICCOUNTS.out.tumor_allelic
  ch_normal_allelic = GATK4_COLLECTALLELICCOUNTS.out.normal_allelic


  GATK4_DENOISEREADCOUNTS(
      ch_counts,
      pon
      )
/*
  GATK4_POLTDENOISEDRATIOS(
      GATK4_DENOISEREADCOUNTS.out.denoised,
      fasta,
      fasta_index,
      fasta_dict
      )
*/
  GATK4_DENOISEREADCOUNTS.out.denoised.concat(
      ch_tumor_allelic,
      ch_normal_allelic
      )
      .groupTuple(by: [0])
      .map {
        meta, v1, v2, v3, v4 ->
        [meta, v1[0], v2, v3, v4, v1[1], v1[2] ]
        }
      .set { ch_seg_input }

  GATK4_MODELSEGMENTS(
      ch_seg_input
      )

  GATK4_DENOISEREADCOUNTS.out.denoised.concat(
      GATK4_MODELSEGMENTS.out.seg
      )
      .groupTuple(by: [0])
      .map{
          meta,v1, v2, v3, v4 ->
          [meta, v1[0], v2, v3, v4, v1[1] ]
      }
      .set { ch_plotseg_input }
/*
  GATK4_POLTMODELSEGMENTS(
      ch_plotseg_input,
      fasta,
      fasta_index,
      fasta_dict
      )
*/
  GATK4_CALLCOPYRATIOSEGMENTS(
      GATK4_MODELSEGMENTS.out.seg
      )

  dataSources = "/data/person/huz/ref/gatk/funcotator/funcotator_dataSources.v1.7.20200521s"
  intervals = "/data/person/huz/DNAseq/UCEC_un_annovar/11_gatk_cnv/ref/agilentv6_wes_S07604514.hg19.interval_list"
/*
  GATK4_FUNCOCATESEGMENTS(
      GATK4_CALLCOPYRATIOSEGMENTS.out.called_seg,
      fasta,
      fasta_index,
      fasta_dict,
      intervals,
      dataSources
      )
*/
  emit:
    seg = GATK4_MODELSEGMENTS.out.seg
    called = GATK4_CALLCOPYRATIOSEGMENTS.out.called_seg

}


workflow {

  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_OV_somatic_bam_channels(it) }
  .set { ch_bam }

/*
ch_input = Channel.fromPath(params.input)
ch_input
.splitCsv( header:true, sep:',' )
.map { create_somatic_bam_channels(it) }
.set { ch_bam }
*/
  genome = params.genome
  fasta = file(params.genomes[genome].fasta)
  fasta_index = file(params.genomes[genome].fasta_fai)
  fasta_dict = file(params.genomes[genome].dict)
  intervals = params.exome ? file(params.genomes[genome].wes_target_bed) : file(params.genomes[genome].wgs_intervals)
  gnomad = file(params.genomes[genome].germline_resource)
  gnomad_idx = file(params.genomes[genome].germline_resource_index)
  variants_for_contamination = file(params.genomes[genome].variants_for_contamination)
  variants_for_contamination_idx = file(params.genomes[genome].variants_for_contamination_index)


  process_intervals = "/data/person/huz/DNAseq/UCEC_un_annovar/11_gatk_cnv/ref/hg19.cnv.preprocessed.interval_list"
  common_sites = "/data/person/huz/DNAseq/UCEC_un_annovar/11_gatk_cnv/ref/common_snps.interval_list"


    hg19_fasta = "/data/share/hg19_wes_ovls/final.fa"
    hg19_fasta_index = "/data/share/hg19_wes_ovls/final.fa.fai"
    hg19_fasta_dict = "/data/share/hg19_wes_ovls/final.dict"
    hg19_process_intervals = "/data/person/huz/DNAseq/OV_long_short/ref/hg19.cnv.preprocessed.interval_list"
    hg19_common_sites = "/data/person/huz/DNAseq/OV_long_short/ref/new.hg19.common_snps.interval_list"

  ch_bam
  .branch {
    meta, bam ->
    single  : bam.size() == 2
    return [ meta, bam[0], bam[1]]
    paired  : bam.size() > 2
    return [ meta, bam[0], bam[1], bam[2], bam[3]]
  }
  .set { ch_paired_bam }

  SOMATIC_CNV(
    ch_paired_bam.paired,
    hg19_fasta,
    hg19_fasta_index,
    hg19_fasta_dict,
    hg19_process_intervals,
    hg19_common_sites

    )
/*
  SOMATIC_CNV(
    ch_paired_bam.paired,
    fasta,
    fasta_index,
    fasta_dict,
    process_intervals
    )
*/
}
