
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process SEQUENZA_UTILS {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_index
    path fasta_dict
    path gc_ref

    output:
    tuple val(meta), path("*.small.seqz.gz"),path("*.small.seqz.gz.tbi"), emit: small_seqz

    script:
    def prefix     = "${meta.id}"

    """
    sequenza-utils bam2seqz \\
    -n ${normal_bam} \\
    -t ${tumor_bam} \\
    --fasta ${fasta} \\
    -gc ${gc_ref} \\
    -o ${prefix}.seqz.gz

    sequenza-utils seqz_binning \\
    --seqz ${prefix}.seqz.gz \\
    -w 50 -o ${prefix}.small.seqz.gz
    """
}

process SEQUENZA_R {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta), path(small_seqz),path(small_seqz_tbi)


    output:
    tuple val(meta), path("Sequenza_Output/*"), emit: sequenza_output

    script:
    def prefix     = "${meta.id}.sequenza"

    """
    #!/usr/bin/env Rscript
    library(sequenza)
    data.file <-  "${small_seqz}"
    extract <- sequenza.extract(data.file, verbose = FALSE)
    CP <- sequenza.fit(extract)
    sequenza.results(sequenza.extract = extract,
                              cp.table = CP,
                              sample.id = "${prefix}",
                              out.dir="Sequenza_Output")
    """
}
