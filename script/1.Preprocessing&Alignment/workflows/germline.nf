nextflow.enable.dsl = 2

include { create_bam_channels } from '../modules/function'

include { GATK4_HAPLOTYPECALLER} from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_GENOMICSDBIIMPORT} from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GATK4_GENOTYPEGVCFS} from '../modules/tools/gatk4' addParams( outdir: params.outdir, publish_enable: true, publish_dir_mode: params.publish_dir_mode)


workflow GERMLINE_CALLING {
    /*
    ch_bam  meta [id, sample, type, ], tumor_bam, tumor_bai, normal_bam, normal_bai
    */
  take:
  ch_bam
  fasta
  fasta_index
  fasta_dict
  intervals
  dbsnp
  dbsnp_index
  use_gvcf

  main:

  GATK4_HAPLOTYPECALLER(
      ch_bam,
      fasta,
      fasta_index,
      fasta_dict,
      intervals,
      dbsnp,
      dbsnp_index,
      use_gvcf
      )
    ch_germ_vcf = GATK4_HAPLOTYPECALLER.out.vcf

    ch_germ_vcf
    .map {
    meta, vcf, tbi ->
    id           = meta.id
    sample       = meta.sample
    sample_type  = meta.sample_type
    patient      = meta.patient
    platform     = meta.platform
    vcf          = "${params.outdir}/${id}.vcf.gz"
    tbi          = "${params.outdir}/${id}.vcf.gz.tbi"
    "${id},${sample},${sample_type},${patient},${platform},${vcf},${tbi}"
    }.set {
    ch_germ_csv
    }

    Channel.from("id,sample,sample_type,patient,platform,vcf,tbi").concat(ch_germ_csv)
    .collectFile(
    name: 'germ_vcf.csv', newLine: true, sort: false, storeDir: "${params.outdir}"
    ).set{ ch_germ_table}

  emit:
    vcf = ch_germ_vcf
    csv = ch_germ_table
  //  DP >=8 VAF >= 0.05 AD >=5 tumor
  //  VAF < 0.01  NAT   strand bias <= 0.95


}


workflow {

  ch_input = Channel.fromPath(params.input)
  ch_input
  .splitCsv( header:true, sep:',' )
  .map { create_germline_bam_channels(it) }
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

  use_gvcf = false

  GERMLINE_CALLING(
    ch_paired_bam.single,
    fasta,
    fasta_index,
    fasta_dict,
    intervals,
    dbsnp,
    dbsnp_index,
    use_gvcf
    )

}
