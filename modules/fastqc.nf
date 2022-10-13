#!/usr/bin/env nextflow

/*
fastQC
*/

process fastqc {
    publishDir params.output, mode: 'copy'

    cpus 16
    time '2h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.fastQC.out/*"

    script:
    """
    zcat ${bam} | head -5000000 > ${bam.baseName}.fastq && \
    mkdir ${bam.baseName}.fastQC.out && \    
    fastqc --threads 16 -o ${bam.baseName}.fastQC.out ${bam.baseName}.fastq
    """
}