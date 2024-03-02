nextflow.enable.dsl = 2

include { create_fastq_channels; create_somatic_bam_channels; create_bam_channels;create_vcf_channels;create_absolute_channels; create_pyclone_channels } from '../modules/function'
include { STRELKA2_SOMATIC } from '../modules/tools/strelka2' addParams( outdir: params.outdir+"/2.manta_strelka2", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { MANTA_SOMATIC } from '../modules/tools/manta' addParams( outdir: params.outdir+"/1.manta", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { STRELKA2VCF2MAF } from '../modules/tools/vcf2maf' addParams( outdir: params.outdir+"/3.variant_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)



workflow STRELKA2 {

  take:
  ch_final_input
  fasta
  fasta_index
  fasta_dict
  intervals
  //gnomad
  //gnomad_idx
  //pon
  //pon_idx
  //variants_for_contamination
  //variants_for_contamination_idx
  //m2_extra_args

  main:

  MANTA_SOMATIC(
    ch_final_input,
    fasta,
    fasta_index,
    fasta_dict,
    intervals
    )

    ch_final_input.concat(
      MANTA_SOMATIC.out.candidatesmallindels_vcf
      )
      .groupTuple(by:[0])
      .map {
        meta, v1, v2, v3, v4 ->
        [meta, v1[0], v2[0], v3, v4, v1[1], v2[1]]
      }
      .set{ ch_input_strelka }

  STRELKA2_SOMATIC(
      ch_input_strelka,
      fasta,
      fasta_index,
      fasta_dict,
      intervals
      )

  emit:
    snv_indels = STRELKA2_SOMATIC.out.snv_indels
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

  vcf2maf_path = "/data/person/huz/my_software/vcf2maf-1.6.21/"
  vep_path = "/data/person/huz/my_software/miniconda3/envs/DNAseq/bin"
  vep_data = "/data/person/huz/ref/vep/homo_sapiens/cache"

  STRELKA2(
    ch_final_input,
    fasta,
    fasta_index,
    fasta_dict,
    intervals
    )

    STRELKA2VCF2MAF(
        STRELKA2.out.snv_indels,
        fasta,
        fasta_index,
        fasta_dict,
        vcf2maf_path,
        vep_path,
        vep_data
    )
}
