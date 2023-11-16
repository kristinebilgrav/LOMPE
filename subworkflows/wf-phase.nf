#!/usr/bin/env nextflow

/*
WORKFLOW forphasing of bam file
*/

// MODULES
include { phase_it ; bamindex } from '../modules/phase'


// WORKFLOW

workflow phasing {
    take : bam_snv_channel

    main: 
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out.phased_bam)

    emit: phased_bam_bai_channel

}