process GROUP_BY_BARCODES {
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.23 bioconda::dnaio=1.2.3 conda-forge::polars=1.38.1 conda-forge::pyarrow=23.0.1 conda-forge::python-duckdb=1.4.4"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:1a408e051b718cb1' :
        'community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:927a9b0310d2951b' }"

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
