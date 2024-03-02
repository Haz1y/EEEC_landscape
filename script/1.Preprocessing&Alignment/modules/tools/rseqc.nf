conda_env = (params.conda_env ? params.conda_env : "bioconda::rseqc=3.0.1")
singularity_url = ""
docker_url = ""

process RSEQC_BAMSTAT {
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
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: txt

    script:
    def prefix   =  "${meta.id}"
    """
    bam_stat.py \\
        -i $bam \\
        > ${prefix}.bam_stat.txt
    """
}

process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: txt

    script:
    def prefix   =  "${meta.id}"
    """
    infer_experiment.py \\
        -i $bam \\
        -r $bed \\
        > ${prefix}.infer_experiment.txt
    """
}

process RSEQC_INNERDISTANCE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*distance.txt"), optional:true, emit: distance
    tuple val(meta), path("*freq.txt")    , optional:true, emit: freq
    tuple val(meta), path("*mean.txt")    , optional:true, emit: mean
    tuple val(meta), path("*.pdf")        , optional:true, emit: pdf
    tuple val(meta), path("*.r")          , optional:true, emit: rscript

    script:
    def prefix   =  "${meta.id}"
        """
        inner_distance.py \\
            -i $bam \\
            -r $bed \\
            -o $prefix \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt

        """
}

process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.junction.bed"), emit: bed
    tuple val(meta), path("*.Interact.bed"), emit: interact_bed
    tuple val(meta), path("*.xls")         , emit: xls
    tuple val(meta), path("*junction.pdf") , emit: pdf
    tuple val(meta), path("*events.pdf")   , emit: events_pdf
    tuple val(meta), path("*.r")           , emit: rscript
    tuple val(meta), path("*.log")         , emit: log

    script:
    def prefix   =  "${meta.id}"
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        2> ${prefix}.junction_annotation.log
    """
}

process RSEQC_JUNCTIONSATURATION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r")  , emit: rscript

    script:
    def prefix   =  "${meta.id}"
    """
    junction_saturation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
    """
}

process RSEQC_READDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.read_distribution.txt"), emit: txt

    script:
    def prefix   =  "${meta.id}"
    """
    read_distribution.py \\
        -i $bam \\
        -r $bed \\
        > ${prefix}.read_distribution.txt
    """
}

process RSEQC_READDUPLICATION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*seq.DupRate.xls"), emit: seq_xls
    tuple val(meta), path("*pos.DupRate.xls"), emit: pos_xls
    tuple val(meta), path("*.pdf")           , emit: pdf
    tuple val(meta), path("*.r")             , emit: rscript

    script:
    def prefix   =  "${meta.id}"
    """
    read_duplication.py \\
        -i $bam \\
        -o $prefix \\
    """
}
