process DORADO_ALIGNER {
    container "nanoporetech/dorado:sha38b4ce849afa13eac8075f0b41cecd30799f169b"
    
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
