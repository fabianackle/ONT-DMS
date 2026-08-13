process REMAP_BAM {
    conda "bioconda::dnaio=1.2.4 bioconda::pysam=0.24.0 bioconda::samtools=1.24"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/dnaio_pysam_samtools:3d162dbb3d2407ff' :
        'community.wave.seqera.io/library/dnaio_pysam_samtools:89063619b55e6ea2' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(alignment), path(barcode_map)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    remap_bam.py \
        --barcode_map ${barcode_map} \
        --input_bam ${alignment} \
        | samtools sort -@ ${task.cpus} --write-index -o ${sample_id}.bam##idx##${sample_id}.bam.bai
    """
}
