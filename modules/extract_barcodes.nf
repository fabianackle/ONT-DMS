process EXTRACT_BARCODES {
    conda "bioconda::cutadapt=5.2"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/cutadapt:5.2--7e96173061760bf9' :
        'community.wave.seqera.io/library/cutadapt:5.2--0bf1fa95d2357cb7' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_barcodes.fasta"), emit: barcodes

    script:
    """
    cutadapt \
        -j ${task.cpus} \
        -g ${params.extract_barcode_adapter} \
        --error-rate 0.1 \
        --revcomp \
        --maximum-length 50 \
        --discard-untrimmed \
        --rename='{header}' \
        --fasta ${fastq_gz} > ${sample_id}_barcodes.fasta
    """
}
