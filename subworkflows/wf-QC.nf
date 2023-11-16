#!/usr/bin/env nextflow

/*
Quality control workflows
*/

include { picard } from './modules/picard'
include { fastqc } from './modules/fastqc'

workflow QC {
    take: 

    main:
    picard(phased_bam_bai_channel)
    fastqc(sample_channel)

    multiQC()
}
