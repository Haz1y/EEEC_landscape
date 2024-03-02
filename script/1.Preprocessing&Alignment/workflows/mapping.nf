#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
input
genomes
genome
exome
outdir
publish_dir_mode
enable_conda
conda_env
*/
include { create_fastq_channels } from '../modules/function'

//include { HEAD } from '../modules/common/head' addParams( outdir: params.outdir, publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { BWA_MEM } from '../modules/tools/bwa' addParams( outdir: params.outdir+"/01_bwa", publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { SAMTOOLS_MERGE } from '../modules/tools/samtools' addParams( outdir: params.outdir+"/01_bwa", publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { GATK4_MARKDUPLICATES} from '../modules/tools/gatk4' addParams( outdir: params.outdir+"/02_markdup", publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { GATK4_BASERECALIBRATOR } from '../modules/tools/gatk4' addParams( outdir: params.outdir+"/03_bqsr", publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { GATK4_APPLYBQSR } from '../modules/tools/gatk4' addParams( outdir: params.outdir+"/03_bqsr", publish_enable: true, publish_dir_mode: params.publish_dir_mode)


workflow MAPPING {

  take:
  ch_fastq
  fasta
  fasta_index
  fasta_dict
  bwa_index
  intervals
  known_sites
  known_sites_tbl
  dbsnp
  dbsnp_index

  main:

  // step 1. align
  //HEAD(ch_fastq)

  //ch_test_fastq = HEAD.out.fastq
  bwa_extra_args = "-K 100000000 -M"
  BWA_MEM(
    ch_fastq,
    bwa_index,
    bwa_extra_args
    )

  BWA_MEM.out.bam
    .map {
      meta, bam ->
      meta.id = meta.sample
      meta.lane = "combined"
      [meta, bam ]
    }
    .groupTuple(by: [0])
    .branch {
      meta, bam ->
      single  : bam.size() == 1
      return [ meta, bam.flatten() ]
      multiple: bam.size() > 1
      return [ meta, bam.flatten() ]
    }.set { ch_mapped_bam }

  SAMTOOLS_MERGE(
    ch_mapped_bam.multiple
    )

  SAMTOOLS_MERGE.out.bam.mix(
    ch_mapped_bam.single
    )
    .set { ch_merged_bam }

  // step 2. mark duplicates
  gatk4_markdup_extra_args = ""

  GATK4_MARKDUPLICATES(
    ch_merged_bam,
    gatk4_markdup_extra_args
    )

  ch_markdup_bam = GATK4_MARKDUPLICATES.out.bam

  // step 3. create recal table

  GATK4_BASERECALIBRATOR(
    ch_markdup_bam,
    fasta,
    fasta_index,
    fasta_dict,
    intervals,
    known_sites,
    known_sites_tbl,
    dbsnp,
    dbsnp_index
    )

  ch_recal_table = ch_markdup_bam.join(GATK4_BASERECALIBRATOR.out.table)

  // STEP 4: RECALIBRATING
  GATK4_APPLYBQSR(
    ch_recal_table,
    fasta,
    fasta_index,
    fasta_dict,
    intervals
    )


  ch_bqsr_bam = GATK4_APPLYBQSR.out.bam

  ch_bqsr_bam
  .map {
    meta, tumor_bam, tumor_bai ->
    id           = meta.id
    sample       = meta.sample
    sample_type  = meta.sample_type
    patient      = meta.patient
    platform     = meta.platform
    bam          = "${params.outdir}/03_bqsr/${id}.bqsr.bam"
    bai          = "${params.outdir}/03_bqsr/${id}.bqsr.bai"
    "${id},${sample},${sample_type},${patient},${platform},${bam},${bai}"
}.set {
    ch_csv
}
Channel.from("id,sample,sample_type,patient,platform,bam,bai").concat(ch_csv)
.collectFile(
    name: 'wgs_bqsr_bam.csv', newLine: true, sort: false, storeDir: "${params.outdir}"
).set{ ch_bqsr_bam_table}

  emit:
     bam = ch_bqsr_bam
     csv = ch_bqsr_bam_table

}


workflow {

  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_fastq_channels(it) }
  .set { ch_fastq }

  genome = params.genome
  bwa_index = file(params.genomes[genome].bwa)
  fasta = file(params.genomes[genome].fasta)
  fasta_index = file(params.genomes[genome].fasta_fai)
  fasta_dict = file(params.genomes[genome].dict)
  intervals = params.exome ? file(params.genomes[genome].wes_intervals) : file(params.genomes[genome].wgs_intervals)
  dbsnp = file(params.genomes[genome].known_indels)
  known_sites = file(params.genomes[genome].known_indels)
  known_sites_tbl = file(params.genomes[genome].known_indels_index)
  dbsnp = file(params.genomes[genome].dbsnp)
  dbsnp_tbl = file(params.genomes[genome].dbsnp_index)

  MAPPING(
    ch_fastq,
    fasta,
    fasta_index,
    fasta_dict,
    bwa_index,
    intervals,
    known_sites,
    known_sites_tbl,
    dbsnp,
    dbsnp_tbl
    )
}
