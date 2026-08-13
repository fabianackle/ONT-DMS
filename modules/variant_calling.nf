process VARIANT_CALLING {
    conda "bioconda::dnaio=1.2.4 conda-forge::biopython=1.88 conda-forge::polars=1.43.2"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/dnaio_biopython_polars:6ad2d6f6836e25da' :
        'community.wave.seqera.io/library/dnaio_biopython_polars:352431a45ac0599c' }"

    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.txt.gz'

    input:
    tuple val(sample_id), path(assembly), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_variants.csv"), emit: variants_db
    tuple val(sample_id), path("${sample_id}_barcodemap.txt.gz"), emit: enrich2_barcodemap

    script:
    """
    variant_calling.py \
        --assembly_path ${assembly} \
        --reference_path ${reference} \
        --sample_id ${sample_id} \
        --barcode_pattern ${params.barcode_5p}\t${params.barcode_3p} \
        --orf_pattern ${params.orf_5p}\t${params.orf_3p} \
        ${params.translate_barcode ? '--translate_barcode' : ''}
        
    gzip ${sample_id}_barcodemap.txt
    """
}
