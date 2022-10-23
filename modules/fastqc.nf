#!/usr/bin/env nextflow

/*
fastQC
*/

process fastqc {
    publishDir params.output, mode: 'copy'

    cpus 16
    time '2h'

    input:
    path(fastq_file)

    output:
    path "${params.sample_id}.fastQC.out/*", emit: QC

    script:
    """
    mkdir ${params.sample_id}.fastQC.out && 
    fastqc --threads 16 -o ${params.sample_id}.fastQC.out ${fastq_file}
    """
}