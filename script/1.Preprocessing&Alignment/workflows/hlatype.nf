nextflow.enable.dsl = 2

include { create_somatic_bam_channels } from '../modules/function'
include { YARA_MAPPING } from '../modules/tools/yara' addParams( outdir: params.outdir +"/01_yara_mapping", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { OPTITYPE } from '../modules/tools/optitype' addParams( outdir: params.outdir +"/02_optitype", publish_enable: true, publish_dir_mode: params.publish_dir_mode)

workflow HLATYPE {

  take:
  ch_bam
  fasta
  fasta_fai
  fasta_dict
  yara_index


  main:
  YARA_MAPPING(
    ch_bam,
    yara_index
    )

  OPTITYPE(
    YARA_MAPPING.out.yara_mapping_bam,
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

    HLATYPE(
      ch_paired_bam.paired,
      fasta,
      fasta_index,
      fasta_dict,
      yara_index
      )
}
