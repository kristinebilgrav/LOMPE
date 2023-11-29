#!/usr/bin/env nextflow

/*
call tandem repeats
*/

process trgt {
    tag "${params.style}:${SampleID}:trgt"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(bam), file(bai)

    output:

    script:
    """

    """

}