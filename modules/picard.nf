#!/usr/bin/env nextflow

/*
picard
*/

process picard {
    tag "${params.style}:${SampleID}:picard"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(SampleID), file(bam), file(bai)

    output:
    tuple val(SampleID), file("${bam.baseName}.wgsmetrics.txt")
    
    shell:
    """
    java -jar \$PICARD_ROOT/picard.jar CollectWgsMetrics --INPUT ${bam} --OUTPUT ${bam.baseName}.wgsmetrics.txt --REFERENCE_SEQUENCE ${params.ref} --COUNT_UNPAIRED true --MINIMUM_BASE_QUALITY 1 --VALIDATION_STRINGENCY SILENT
    """
}