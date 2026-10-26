process GROUP_BY_BARCODES {
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.24 bioconda::dnaio=1.2.4 conda-forge::polars=2.1.0 conda-forge::pyarrow=25.0.0 conda-forge::python-duckdb=2.0.0"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/bwa_samtools_dnaio_polars_pruned:0a127e00d4efd86e' :
        'community.wave.seqera.io/library/bwa_samtools_dnaio_polars_pruned:9ed0d29bd48be8c0' }"

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
