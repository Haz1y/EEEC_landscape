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

include { FASTP } from '../modules/tools/fastp' addParams( outdir: params.outdir, publish_enable: params.publish_enable, publish_dir_mode: params.publish_dir_mode)


workflow FASTQCLEAN {

    take:
    ch_fastq

    main:

    fastp_extra_args = "--length_required 20"
    FASTP(
        ch_fastq,
        fastp_extra_args
        )

    FASTP.out.reads
    .map {
        meta, reads ->
        id           = meta.id
        sample       = meta.sample
        sample_type  = meta.sample_type
        patient      = meta.patient
        platform     = meta.platform
        lane         = meta.lane
        fastq1 = "${params.outdir}/${id}.clean.R1.fastq.gz"
        fastq2 = "${params.outdir}/${id}.clean.R2.fastq.gz"
        "${id},${sample},${sample_type},${patient},${platform},${lane},${fastq1},${fastq2}"
        }
    .set {
        ch_csv
    }

    Channel.from("id,sample,sample_type,patient,platform,lane,forward_file_name,reverse_file_name").concat(ch_csv)
    .collectFile(
        name: 'clean_fastq.csv', newLine: true, sort: false, storeDir: "${params.outdir}"
    )
    
    emit:
    reads = FASTP.out.reads
    csv = "${params.outdir}/clean_fastq.csv"

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
  known_sites = file(params.genomes[genome].known_indels)
  known_sites_tbl = file(params.genomes[genome].known_indels_index)
  dbsnp = file(params.genomes[genome].dbsnp)
  dbsnp_tbl = file(params.genomes[genome].dbsnp_index)

  FASTQCLEAN(
    ch_fastq
    )
}
