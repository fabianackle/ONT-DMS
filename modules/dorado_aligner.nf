process DORADO_ALIGNER {
    container "nanoporetech/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e"
    
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
