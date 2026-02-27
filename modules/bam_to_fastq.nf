process BAM_TO_FASTQ {
    conda "bioconda::samtools=1.23 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/samtools_pigz:8e0bfad3948505f5' :
        'community.wave.seqera.io/library/samtools_pigz:a80e8586599cb352' }"

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
