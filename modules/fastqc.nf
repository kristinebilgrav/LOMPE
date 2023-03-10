#!/usr/bin/env nextflow

/*
fastQC
*/

process fastqc {
    tag "${params.style}:${SampleID}:fastQC"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(fastq_file)

    output:
    tuple val(SampleID), file("${SampleID}.fastQC.out/*"), emit: QC

    script:
    """
    mkdir ${SampleID}.fastQC.out && 
    fastqc --threads 16 -o ${SampleID}.fastQC.out ${fastq_file}
    """
}