nextflow.enable.dsl = 2

include { create_somatic_bam_channels } from '../modules/function'
include { SEQUENZA_UTILS } from '../modules/tools/sequenza' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { SEQUENZA_R } from '../modules/tools/sequenza' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)

workflow SEQUENZA {

  take:
  ch_bam
  fasta
  fai
  dict
  gc_ref


  main:
  SEQUENZA_UTILS(
    ch_bam,
    fasta,
    fai,
    dict,
    gc_ref
    )

  SEQUENZA_R(
    SEQUENZA_UTILS.out.small_seqz
    )

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
    gc_ref = "/data/person/huz/DNAseq/Sequenza_UCEC/hg19.gc50Base.wig.gz"

    ch_bam
    .branch {
      meta, bam ->
      single  : bam.size() == 2
      return [ meta, bam[0], bam[1]]
      paired  : bam.size() > 2
      return [ meta, bam[0], bam[1], bam[2], bam[3]]
    }
    .set { ch_paired_bam }

    SEQUENZA(
      ch_paired_bam.paired,
      fasta,
      fasta_index,
      fasta_dict,
      gc_ref
      )
}
