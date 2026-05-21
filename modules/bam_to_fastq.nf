process BAM_TO_FASTQ {
    conda "bioconda::samtools=1.23.1 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/samtools_pigz:f21d48b24bf42a88' :
        'community.wave.seqera.io/library/samtools_pigz:81f4c8cda5824ba2' }"

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
