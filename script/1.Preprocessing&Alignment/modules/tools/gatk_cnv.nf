conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process GATK4_COLLECTREADCOUNTS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_index
    path fasta_dict
    path process_intervals

    output:
    tuple val(meta), path("tumor/${meta.id}.tumor.clean_counts.hdf5"), emit: tumor_hdf5
    tuple val(meta), path("normal/${meta.id}.normal.clean_counts.hdf5"), emit: normal_hdf5

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    mkdir tumor
    mkdir normal
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CollectReadCounts \\
        -I ${tumor_bam} \\
        -L ${process_intervals} \\
        -R ${fasta} \\
        --format HDF5 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output tumor/${prefix}.tumor.clean_counts.hdf5

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CollectReadCounts \\
        -I ${normal_bam} \\
        -L ${process_intervals} \\
        -R ${fasta} \\
        --format HDF5 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output normal/${prefix}.normal.clean_counts.hdf5
    """
}

process GATK4_CREATEREADCOUNTPANELOFNORMALS {
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    path(normal_counts)

    output:
    path("cnvponM.pon.hdf5"), emit: pon

    script:
    pon_inputs = normal_counts.collect({v -> "--input $v "}).join(" ")
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CreateReadCountPanelOfNormals \\
        --minimum-interval-median-percentile 5.0 \\
        --output cnvponM.pon.hdf5 \\
        ${pon_inputs}
    """
}

process GATK4_COLLECTALLELICCOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_index
    path fasta_dict
    path common_sites

    output:
    tuple val(meta), path("${meta.id}.tumor.allelicCounts.tsv"), emit: tumor_allelic
    tuple val(meta), path("${meta.id}.normal.allelicCounts.tsv"), emit: normal_allelic

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CollectAllelicCounts \\
        -L ${common_sites} \\
        --input ${tumor_bam} \\
        --reference ${fasta} \\
        --minimum-base-quality 20 \\
        --output ${prefix}.tumor.allelicCounts.tsv

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CollectAllelicCounts \\
        -L ${common_sites} \\
        --input ${normal_bam} \\
        --reference ${fasta} \\
        --minimum-base-quality 20 \\
        --output ${prefix}.normal.allelicCounts.tsv
    """
}

process GATK4_DENOISEREADCOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_counts), path(normal_counts)
    path pon

    output:
    tuple val(meta), path("${meta.id}.tumor.standardizedCR.tsv"), path("${meta.id}.tumor.denoisedCR.tsv"), path("${meta.id}.normal.standardizedCR.tsv"), path("${meta.id}.normal.denoisedCR.tsv"), emit: denoised

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        DenoiseReadCounts \\
        -I ${tumor_counts} \\
        --count-panel-of-normals ${pon} \\
        --standardized-copy-ratios ${prefix}.tumor.standardizedCR.tsv \\
        --denoised-copy-ratios ${prefix}.tumor.denoisedCR.tsv \

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        DenoiseReadCounts \\
        -I ${normal_counts} \\
        --count-panel-of-normals ${pon} \\
        --standardized-copy-ratios ${prefix}.normal.standardizedCR.tsv \\
        --denoised-copy-ratios ${prefix}.normal.denoisedCR.tsv \
    """
}


process GATK4_MODELSEGMENTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_standard), path(tumor_denoise), path(normal_standard), path(normal_denoise), path(tumor_allelic), path(normal_allelic)

    output:
    tuple val(meta), path("segments/*"), emit: seg

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        ModelSegments \\
        --denoised-copy-ratios ${tumor_denoise} \\
        --allelic-counts ${tumor_allelic} \\
        --normal-allelic-counts ${normal_allelic} \\
        --output segments \\
        --output-prefix ${prefix}.tumor

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        ModelSegments \\
        --denoised-copy-ratios ${normal_denoise} \\
        --allelic-counts ${normal_allelic} \\
        --output segments \\
        --output-prefix ${prefix}.normal
    """
}

process GATK4_CALLCOPYRATIOSEGMENTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path("segments/*")

    output:
    tuple val(meta), path("segments/${meta.id}.tumor.called.seg"),path("segments/${meta.id}.normal.called.seg"), emit: called_seg

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CallCopyRatioSegments \\
        -I segments/${prefix}.tumor.cr.seg \\
        -O segments/${prefix}.tumor.called.seg

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CallCopyRatioSegments \\
        -I segments/${prefix}.normal.cr.seg \\
        -O segments/${prefix}.normal.called.seg
    """
}

process GATK4_POLTDENOISEDRATIOS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_standard), path(tumor_denoise), path(normal_standard), path(normal_denoise)
    path(fasta)
    path(fasta_index)
    path(fasta_dict)

    output:
    tuple val(meta), path("plots_denoise/*"), emit: plots_denoise

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        PlotDenoisedCopyRatios \\
        --standardized-copy-ratios ${tumor_standard} \\
        --denoised-copy-ratios ${tumor_denoise} \\
        --sequence-dictionary ${fasta_dict} \\
        --output plots_denoise \\
        --output-prefix ${prefix}.tumor

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        PlotDenoisedCopyRatios \\
        --standardized-copy-ratios ${normal_standard} \\
        --denoised-copy-ratios ${normal_denoise} \\
        --sequence-dictionary ${fasta_dict} \\
        --output plots_denoise \\
        --output-prefix ${prefix}.normal
    """
}

process GATK4_POLTMODELSEGMENTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_standard), path(tumor_denoise), path(normal_standard), path(normal_denoise),  path("segments/*")
    path(fasta)
    path(fasta_index)
    path(fasta_dict)

    output:
    tuple val(meta), path("plots_model_seg/*"), emit: plots_model_seg

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        PlotModeledSegments \\
        --denoised-copy-ratios ${tumor_denoise} \\
        --allelic-counts segments/${prefix}.tumor.hets.tsv \\
        --segments segments/${prefix}.tumor.modelFinal.seg \\
        --sequence-dictionary ${fasta_dict} \\
        --output plots_model_seg \\
        --output-prefix ${prefix}.tumor

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        PlotModeledSegments \\
        --denoised-copy-ratios ${normal_denoise} \\
        --allelic-counts segments/${prefix}.normal.hets.tsv \\
        --segments segments/${prefix}.normal.modelFinal.seg \\
        --sequence-dictionary ${fasta_dict} \\
        --output plots_model_seg \\
        --output-prefix ${prefix}.normal
    """
}

process GATK4_FUNCOCATESEGMENTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
      mode: params.publish_dir_mode,
      enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/DNAseq" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_called_seg),path(normal_called_seg)
    path(fasta)
    path(fasta_index)
    path(fasta_dict)
    path(intervals)
    val (dataSources)

    output:
    tuple val(meta), path("funcotated/*"), emit: funcotated_seg

    script:
    def prefix    = "${meta.id}"
    // def normal_command = meta.normal ? "-matched ${normal_pileup}" : ""
    """
    mkdir funcotated

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        FuncotateSegments \\
        --data-sources-path ${dataSources} \\
        --ref-version hg19 \\
        --output-file-format SEG \\
        -R ${fasta} \\
        --segments ${tumor_called_seg} \\
        -O funcotated/${prefix}.tumor.funcotated.seg \\
        --transcript-list funcotated/${prefix}.tumor.tx_list.txt

    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        FuncotateSegments \\
        --data-sources-path ${dataSources} \\
        --ref-version hg19 \\
        --output-file-format SEG \\
        -R ${fasta} \\
        --segments ${normal_called_seg} \\
        -O funcotated/${prefix}.normal.funcotated.seg \\
        --transcript-list funcotated/${prefix}.normal.tx_list.txt
    """
}
