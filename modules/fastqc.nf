#!/usr/bin/env nextflow

/*
fastQC
*/

process fastqc {
    publishDir params.outdir, mode: 'copy'

    cpus 16
    time '2h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.QC.out"

    script:
    """
    zcat ${bam} | head -5000000 > ${bam.baseName}.fastq
    mkdir ${bam.baseName}
    fastqc --threads 16 -o ${bam.baseName}.fastQC.out ${bam.baseName}.fastq
    """
}