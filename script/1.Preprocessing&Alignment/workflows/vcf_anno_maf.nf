nextflow.enable.dsl = 2

include { create_vcf_channels } from '../modules/function'

include { VCF2MAF} from '../modules/tools/vcf2maf' addParams( outdir: params.outdir+"/09_variant_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)




workflow {

  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_vcf_channels(it) }
  .set { ch_vcf }

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

  vcf2maf_path = "/data/person/huz/my_software/vcf2maf-1.6.21/"
  vep_path = "/data/person/huz/my_software/miniconda3/envs/DNAseq/bin"
  vep_data = "/data/person/huz/ref/vep/homo_sapiens/cache"

  VCF2MAF(
    ch_vcf,
    fasta,
    fasta_index,
    fasta_dict,
    vcf2maf_path,
    vep_path,
    vep_data
    )

}
