conda_env = (params.conda_env ? params.conda_env : "bioconda::rsem=1.3.3 bioconda::star=2.7.6a bioconda::samtools=1.10")
singularity_url = ""
docker_url = ""

process RSEM_CALCULATEEXPRESSION {
    tag "$meta.id"
    label 'process_high'
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
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat

    script:
    prefix         ="${meta.id}.rsem"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--strandedness forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--strandedness reverse'
    }
    def paired_end = meta.single_end ? "" : "--paired-end"
    """
    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        $strandedness \\
        $paired_end \\
        -no-bam-output \\
        --alignments \\
        --estimate-rspd \\
        --seed 1 \\
        $reads \\
        $index/RSEMIndex \\
        $prefix
    """
}


process RSEM_MERGE_COUNTS {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path ('genes/*')
    path ('isoforms/*')

    output:
    path "rsem.merged.gene_counts.tsv"      , emit: counts_gene
    path "rsem.merged.gene_tpm.tsv"         , emit: tpm_gene
    path "rsem.merged.transcript_counts.tsv", emit: counts_transcript
    path "rsem.merged.transcript_tpm.tsv"   , emit: tpm_transcript

    script:
    """
    mkdir -p tmp/genes
    cut -f 1,2 `ls ./genes/* | head -n 1` > gene_ids.txt
    for fileid in `ls ./genes/*`; do
        samplename=`basename \$fileid | sed s/\\.genes.results\$//g`
        echo \$samplename > tmp/genes/\${samplename}.counts.txt
        cut -f 5 \${fileid} | tail -n+2 >> tmp/genes/\${samplename}.counts.txt
        echo \$samplename > tmp/genes/\${samplename}.tpm.txt
        cut -f 6 \${fileid} | tail -n+2 >> tmp/genes/\${samplename}.tpm.txt
    done

    mkdir -p tmp/isoforms
    cut -f 1,2 `ls ./isoforms/* | head -n 1` > transcript_ids.txt
    for fileid in `ls ./isoforms/*`; do
        samplename=`basename \$fileid | sed s/\\.isoforms.results\$//g`
        echo \$samplename > tmp/isoforms/\${samplename}.counts.txt
        cut -f 5 \${fileid} | tail -n+2 >> tmp/isoforms/\${samplename}.counts.txt
        echo \$samplename > tmp/isoforms/\${samplename}.tpm.txt
        cut -f 6 \${fileid} | tail -n+2 >> tmp/isoforms/\${samplename}.tpm.txt
    done

    paste gene_ids.txt tmp/genes/*.counts.txt > rsem.merged.gene_counts.tsv
    paste gene_ids.txt tmp/genes/*.tpm.txt > rsem.merged.gene_tpm.tsv
    paste transcript_ids.txt tmp/isoforms/*.counts.txt > rsem.merged.transcript_counts.tsv
    paste transcript_ids.txt tmp/isoforms/*.tpm.txt > rsem.merged.transcript_tpm.tsv
    """
}

process RSEM_PREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? "${conda_env}" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    }

    input:
    path fasta
    path gtf

    output:
    path "*transcripts.fa", emit: transcript_fasta

    script:
        """
        rsem-prepare-reference \\
            --num-threads $task.cpus \\
            --gtf $gtf \\
            $fasta \\
            RSEMIndex
        """
}
