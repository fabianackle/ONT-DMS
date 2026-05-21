process CONSEQUENCE_CALLING {
    conda "bioconda::dnaio=1.2.3 bioconda::bcftools=1.23.1"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/bcftools_dnaio:cc9d8b60c21ae82f' :
        'community.wave.seqera.io/library/bcftools_dnaio:aa359363d1c71dda' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(variants), path(references)

    output:
    tuple val(sample_id), path("${sample_id}_consequences.csv"), emit: consequences

    script:
    """
    write_gff.py \
        --references ${references} \
        --orf_5p ${params.orf_5p} \
        --orf_3p ${params.orf_3p} \
        --output ${sample_id}.gff

    bcftools csq \
        -f ${references} \
        -g ${sample_id}.gff \
        ${variants} \
        > ${sample_id}_consequences.vcf

    bcftools query \
        -f '[%CHROM\t%QUAL\t%FILTER\t%TBCSQ{0}\n]' \
        ${sample_id}_consequences.vcf \
        | parse_consequences.py \
        > ${sample_id}_consequences.csv
    """
}
