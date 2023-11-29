#!/usr/bin/env nextflow

/*
WORKFLOW for SNV calling and filtering
*/

// MODULES

include { sniff } from '../modules/sniffles'
include { run_vep } from '../modules/annotate'

workflow sv_calling {
    take: 
    phased_bam_bai_channel

    main:
    sniff(phased_bam_bai_channel)
}
