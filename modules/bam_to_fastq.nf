process BAM_TO_FASTQ {
    conda "bioconda::samtools=1.25 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/samtools_pigz:f09f6e6f5f106efb' :
        'community.wave.seqera.io/library/samtools_pigz:d61990c90b92de30' }"

    tag "${basecalled.baseName}"

    input:
    tuple val(sample_id), path(basecalled)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq_gz

    script:
    """
    samtools fastq \
        --threads ${task.cpus} \
        ${basecalled} \
        | pigz -p ${task.cpus} > ${sample_id}.fastq.gz
    """
}
