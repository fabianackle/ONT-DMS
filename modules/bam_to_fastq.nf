process BAM_TO_FASTQ {
    conda "bioconda::samtools=1.24 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/samtools_pigz:9d409f843d9e96ff' :
        'community.wave.seqera.io/library/samtools_pigz:261c0b94a7fc9366' }"

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
