#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { create_fastq_channels; create_maf_channels; create_somatic_bam_channels; create_bam_channels;create_vcf_channels;create_absolute_channels; create_pyclone_channels } from './modules/function'

include { FASTQCLEAN } from './workflows/fastq_clean' addParams( outdir: params.outdir+"/01_clean", publish_enable: true, publish_report_only: true, publish_dir_mode: params.publish_dir_mode)
include { MAPPING } from './workflows/mapping' addParams( outdir: params.outdir+"/02_mapping", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { CREATE_SOMATIC_BAM_TABLE } from './workflows/python' addParams( outdir: params.outdir+"/02_mapping", publish_enable: false, publish_dir_mode: params.publish_dir_mode)
include { MUTECT2_PAIRED } from './workflows/mutect2' addParams( outdir: params.outdir+"/03_variant_calling", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { GERMLINE_CALLING } from './workflows/germline' addParams( outdir: params.outdir+"/05_germ_variant_calling", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { FILTER_MUTEC2_VCF } from './workflows/python' addParams( outdir: params.outdir+"/03_variant_calling", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { ANNOVAR } from './modules/tools/annovar' addParams( outdir: params.outdir+"/04_variant_annotate", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { CNVKIT } from './modules/tools/cnvkit' addParams( outdir: params.outdir+"/06_somatic_cnv_calling", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { SEQUENZA } from './workflows/sequenza' addParams( outdir: params.outdir+"/07_purity_sequenza", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { CNV_FACETS } from './modules/tools/facets' addParams( outdir: params.outdir+"/08_facets", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { VCF2MAF } from './modules/tools/vcf2maf' addParams( outdir: params.outdir+"/09_variant_maf", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { DOABSOLUTE } from './modules/tools/absolute' addParams( outdir: params.outdir+"/12_gatk_cnv_absolute", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { SOMATIC_CNV } from './workflows/somatic_cnv' addParams( outdir: params.outdir+"/11_gatk_cnv", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { PYCLONE } from './modules/tools/pyclone' addParams( outdir: params.outdir+"/13_pyclone", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { HLATYPE } from './workflows/hlatype' addParams( outdir: params.outdir+"/14_hlatype", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { COPYWRITER } from './modules/tools/copywriteR' addParams( outdir: params.outdir+"/15_Copywrite", publish_enable: true, publish_dir_mode: params.publish_dir_mode)
include { PURBAYES } from './modules/tools/purbayes' addParams( outdir: params.outdir+"/16_Purbayes", publish_enable: true, publish_dir_mode: params.publish_dir_mode)


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
    if (mode == "mapping") {
        mapping = true
        somatic = false
        germline = false
        annotate = false
        cnv = false
        sequenza =false
        vcf2maf = false
    } else if (mode == "somatic") {
        mapping = true
        somatic = true
        germline = false
        annotate = true
        cnv = false
        sequenza =false
        vcf2maf =false
    } else if (mode == "all") {
        mapping = true
        somatic = true
        germline = true
        annotate = true
        cnv = false
        sequenza =false
        vcf2maf = false
    } 

    if (mapping) {
        ch_input = Channel.fromPath(params.input)

        ch_input
        .splitCsv( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { ch_fastq }

        FASTQCLEAN(
            ch_fastq
        )

        MAPPING(
            FASTQCLEAN.out.reads,
            fasta,
            fasta_index,
            fasta_dict,
            bwa_index,
            intervals,
            known_sites,
            known_sites_tbl,
            dbsnp,
            dbsnp_index
        )
    }

    if (somatic) {
        CREATE_SOMATIC_BAM_TABLE(
            MAPPING.out.csv
            )

        CREATE_SOMATIC_BAM_TABLE.out.table
        .splitCsv( header:true, sep:',' )
        .map { create_somatic_bam_channels(it) }
        .set { ch_bam }

        ch_bam
        .branch {
          meta, bam ->
          single  : bam.size() == 2
          return [ meta, bam[0], bam[1]]
          paired  : bam.size() > 2
          return [ meta, bam[0], bam[1], bam[2], bam[3]]
        }
        .set { ch_paired_bam }

        MUTECT2_PAIRED(
            ch_paired_bam.paired,
            fasta,
            fasta_index,
            fasta_dict,
            intervals,
            gnomad,
            gnomad_idx,
            variants_for_contamination,
            variants_for_contamination_idx,
            m2_extra_args
        )

        min_tumor_dp = 8
        min_normal_dp = 8
        min_tumor_ad = 5
        min_tumor_af = 0.05
        max_normal_af = 0.02

        FILTER_MUTEC2_VCF(
            MUTECT2_PAIRED.out.vcf,
            min_tumor_dp,
            min_normal_dp,
            min_tumor_ad,
            min_tumor_af,
            max_normal_af
        )
    }

    if (germline) {

      MAPPING.out.csv
      .splitCsv( header:true, sep:',' )
      .map { create_bam_channels(it) }
      .set { ch_germline_bam }

      use_gvcf = false

      GERMLINE_CALLING(
          ch_germline_bam,
          fasta,
          fasta_index,
          fasta_dict,
          intervals,
          dbsnp,
          dbsnp_index,
          use_gvcf
      )
    }

    if (annotate) {
        annovar_path = "/data/software/annovar/"

        ANNOVAR(
            FILTER_MUTEC2_VCF.out.vcf,
            fasta,
            fasta_index,
            fasta_dict,
            genome,
            annovar_path
        )
    }
    if (vcf2maf) {
        vcf2maf_path = "/data/person/huz/my_software/vcf2maf-1.6.21/"
        vep_path = "/data/person/huz/my_software/miniconda3/envs/DNAseq/bin"
        vep_data = "/data/person/huz/ref/vep/homo_sapiens/cache"
        ch_input = Channel.fromPath(params.input)

        ch_input
        .splitCsv( header:true, sep:',' )
        .map { create_vcf_channels(it) }
        .set { ch_vcf }

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

    if (cnv) {
        refFlat = " "

        ch_input = Channel.fromPath(params.input)

        ch_input
        .splitCsv( header:true, sep:',' )
        .map { create_somatic_bam_channels(it) }
        .set { ch_bam }

        ch_bam
        .branch {
          meta, bam ->
          single  : bam.size() == 2
          return [ meta, bam[0], bam[1]]
          paired  : bam.size() > 2
          return [ meta, bam[0], bam[1], bam[2], bam[3]]
        }
        .set { ch_paired_bam }

        CNVKIT(
            ch_paired_bam.paired,
            fasta,
            fasta_index,
            fasta_dict,
            intervals,
            refFlat,
            params.exome
        )
    }
    if (absolute) {
        ch_input = Channel.fromPath(params.input)

        ch_input
        .splitCsv( header:true, sep:',' )
        .map { create_absolute_channels(it) }
        .set { ch_absolute }


        DOABSOLUTE(
            ch_absolute
        )
    }
