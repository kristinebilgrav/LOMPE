#!/usr/bin/env nextflow

/*
picard
*/

process picard {
    publishDir params.output, mode:'copy'
    cpus 2
    time '2h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.wgsmetrics.txt"
    
    shell:
    """
    java -jar \$PICARD_ROOT/picard.jar CollectWgsMetrics --INPUT ${bam} --OUTPUT ${bam.baseName}.wgsmetrics.txt --REFERENCE_SEQUENCE ${params.ref} --COUNT_UNPAIRED true --MINIMUM_BASE_QUALITY 1 --VALIDATION_STRINGENCY SILENT
    """
}