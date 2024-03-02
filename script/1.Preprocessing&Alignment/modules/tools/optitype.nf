
conda_env = (params.conda_env? params.conda_env : "bioconda::bcftools perl-getopt-long perl-pod-usage perl-file-spec")
singularity_url = ""
docker_url = ""

process OPTITYPE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        enabled: params.publish_enable

    conda (params.enable_conda ? "/data/person/huz/my_software/miniconda3/envs/pyclone" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${singularity_url}"
    } else {
        container "${docker_url}"
    }

    input:
    tuple val(meta), path(tumor_mapped_bam), path(tumor_mapped_bai), path(normal_mapped_bam), path(normal_mapped_bai)

    output:
    tuple val(meta), path("output/*"), emit: optitype_output

    script:
    def prefix     = "${meta.patient}"
    def id     = "${meta.id}"

    """
    # Create a config for OptiType on a per sample basis with task.ext.args2


    # Run the actual OptiType typing with args
    OptiTypePipeline.py -i ${tumor_mapped_bam}  -e 1 -b 0.009 \\
        -p ${id} --dna --outdir output

    OptiTypePipeline.py -i ${normal_mapped_bam}  -e 1 -b 0.009 \\
        -p ${prefix}_N --dna --outdir output
    """
}
