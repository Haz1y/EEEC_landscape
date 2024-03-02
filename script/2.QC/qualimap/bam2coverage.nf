#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { create_fastq_channels; create_maf_channels; create_somatic_bam_channels; create_bam_channels;create_vcf_channels;create_absolute_channels; create_pyclone_channels } from './modules/function'


include { QUALIMAP_BAMQC } from './modules/tools/qualimap' addParams( outdir: params.outdir+"/04_qualimap_output", publish_enable: true, publish_dir_mode: params.publish_dir_mode)


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


      ch_input = Channel.fromPath(params.input)

      ch_input
      .splitCsv( header:true, sep:',' )
      .map { create_bam_channels(it) }
      .set { ch_bam }
      QUALIMAP_BAMQC(
        ch_bam,
        intervals
        )

}
