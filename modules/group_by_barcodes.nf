process GROUP_BY_BARCODES {
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.23.1 bioconda::dnaio=1.2.3 conda-forge::polars=1.40.1 conda-forge::pyarrow=24.0.0 conda-forge::python-duckdb=1.5.3"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:af08e02186cb2009' :
        'community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:8056255e5d60d4e1' }"

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
