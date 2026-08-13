process DORADO_ALIGNER {
    container "nanoporetech/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(basecalled), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: alignment

    script:
    """
    dorado aligner ${reference} ${basecalled} > ${sample_id}.aligned.bam
    """
}
