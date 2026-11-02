process VARIANT_CALLING {
    conda "bioconda::dnaio=1.2.4 conda-forge::biopython=1.88 conda-forge::gzip=1.14 conda-forge::polars=2.1.0"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/dnaio_biopython_gzip_polars:598924ba7ca09670' :
        'community.wave.seqera.io/library/dnaio_biopython_gzip_polars:f7a5574515f59aa6' }"

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
