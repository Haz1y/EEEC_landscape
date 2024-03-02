
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process YARA_INDEX {
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
    path fasta
    path fasta_index
    path fasta_dict

    output:
    path "yara"        , emit: yara_index
    path "versions.yml", emit: versions

    script:

    """
    mkdir yara

    yara_indexer \\
        $fasta \\
        -o "yara"

    mv *.{lf,rid,sa,txt}.* yara

    cp $fasta yara/yara.fasta

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        yara: \$(echo \$(yara_indexer --version 2>&1) | sed 's/^.*yara_indexer version: //; s/ .*\$//')
    END_VERSIONS

    """
}

process YARA_MAPPING {
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
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path yara_index

    output:
    tuple val(meta), path("${meta.patient}.tumor_mapped_sort_{1,2}.bam"),path("${meta.patient}.tumor_mapped_sort_{1,2}.bam.bai"),path("${meta.patient}.normal_mapped_sort_{1,2}.bam"),path("${meta.patient}.normal_mapped_sort_{1,2}.bam.bai"), emit: yara_mapping_bam

    script:
    def prefix     = "${meta.patient}"

    """
    samtools view -@ ${task.cpus} -bh -f 0x40 ${tumor_bam} > tumor_output_1.bam
    samtools view -@ ${task.cpus} -bh -f 0x80 ${tumor_bam} > tumor_output_2.bam
    samtools bam2fq -@ ${task.cpus} tumor_output_1.bam > tumor_output_1.fastq
    samtools bam2fq -@ ${task.cpus} tumor_output_2.bam > tumor_output_2.fastq
    yara_mapper -e 3 -t ${task.cpus} -f bam ${yara_index}/hla_reference_dna tumor_output_1.fastq tumor_output_2.fastq > tumor_output.bam
    samtools view -@ ${task.cpus} -bh -F 4 -f 0x40 -b1 tumor_output.bam > ${prefix}.tumor_mapped_1.bam
    samtools view -@ ${task.cpus} -bh -F 4 -f 0x80 -b1 tumor_output.bam > ${prefix}.tumor_mapped_2.bam
    samtools sort -@ ${task.cpus} -m 4G ${prefix}.tumor_mapped_1.bam > ${prefix}.tumor_mapped_sort_1.bam
    samtools sort -@ ${task.cpus} -m 4G ${prefix}.tumor_mapped_2.bam > ${prefix}.tumor_mapped_sort_2.bam
    samtools index -@ ${task.cpus} ${prefix}.tumor_mapped_sort_1.bam
    samtools index -@ ${task.cpus} ${prefix}.tumor_mapped_sort_2.bam

    samtools view -@ ${task.cpus} -bh -f 0x40 ${normal_bam} > normal_output_1.bam
    samtools view -@ ${task.cpus} -bh -f 0x80 ${normal_bam} > normal_output_2.bam
    samtools bam2fq -@ ${task.cpus} normal_output_1.bam > normal_output_1.fastq
    samtools bam2fq -@ ${task.cpus} normal_output_2.bam > normal_output_2.fastq
    yara_mapper -e 3 -t ${task.cpus} -f bam ${yara_index}/hla_reference_dna normal_output_1.fastq normal_output_2.fastq > normal_output.bam
    samtools view -@ ${task.cpus} -bh -F 4 -f 0x40 -b1 normal_output.bam > ${prefix}.normal_mapped_1.bam
    samtools view -@ ${task.cpus} -bh -F 4 -f 0x80 -b1 normal_output.bam > ${prefix}.normal_mapped_2.bam
    samtools sort -@ ${task.cpus} -m 4G ${prefix}.normal_mapped_1.bam > ${prefix}.normal_mapped_sort_1.bam
    samtools sort -@ ${task.cpus} -m 4G ${prefix}.normal_mapped_2.bam > ${prefix}.normal_mapped_sort_2.bam
    samtools index -@ ${task.cpus} ${prefix}.normal_mapped_sort_1.bam
    samtools index -@ ${task.cpus} ${prefix}.normal_mapped_sort_2.bam

    rm -rf tumor_output*
    rm -rf normal_output*

    """
}
