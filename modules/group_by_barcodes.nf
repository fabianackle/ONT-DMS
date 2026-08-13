process GROUP_BY_BARCODES {
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.24 bioconda::dnaio=1.2.4 conda-forge::polars=1.43.2 conda-forge::pyarrow=25.0.0 conda-forge::python-duckdb=1.5.5"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:82e47043e312ada4' :
        'community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:8841e9b2d39b827e' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(barcodes), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_references.fasta"), emit: references
    tuple val(sample_id), path("${sample_id}_mapped_reads_filtered.csv"), emit: barcode_map
    tuple val(sample_id), path("${sample_id}_reads.csv"), path("${sample_id}_hq_barcodes.csv"), path("${sample_id}_mapped_reads.csv"), path("${sample_id}_mapped_reads_filtered.csv"), emit: csv

    script:
    """
    group_by_barcodes.py \
        --sample_id ${sample_id} \
        --reference_seq ${reference} \
        --barcodes ${barcodes} \
        --barcode_regex "${params.barcode_regex}" \
        --barcode_min_coverage ${params.barcode_min_coverage} \
        --threads ${task.cpus}
    """
}
