nextflow.enable.dsl = 2

include { create_germ_vcf_channels } from '../modules/function'

include { GATK4_TOTAL_GERMLINE_VQSR} from '../modules/tools/gatk4' addParams( outdir: params.outdir+"/26_vqsr_germ_vcf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GERMLINE_ANNOVAR} from '../modules/tools/annovar' addParams( outdir: params.outdir+"/27_annovar_germ_vcf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)




workflow {

  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_germ_vcf_channels(it) }
  .set { ch_germ_vcf }

  genome = params.genome
  fasta = file(params.genomes[genome].fasta)
  fasta_index = file(params.genomes[genome].fasta_fai)
  fasta_dict = file(params.genomes[genome].dict)
  intervals = params.exome ? file(params.genomes[genome].wes_target_bed) : file(params.genomes[genome].wgs_intervals)

  known_indels = file(params.genomes[genome].known_indels)
  known_indels_index = file(params.genomes[genome].known_indels_index)


  GATK4_TOTAL_GERMLINE_VQSR(
    ch_germ_vcf,
    fasta,
    fasta_index,
    fasta_dict,
    known_indels,
    known_indels_index
    )

  annovar_path = "/data/software/annovar/"

  GERMLINE_ANNOVAR(
    GATK4_TOTAL_GERMLINE_VQSR.out.vqsr_vcf,
    fasta,
    fasta_index,
    fasta_dict,
    annovar_path
    )
}