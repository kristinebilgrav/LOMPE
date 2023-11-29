#!/usr/bin/env nextflow

/*
WORKFLOW start at fastq files
*/

// MODULES

if (params.style == 'ont' && params.file == 'fastq') {
    include {cat} from '../modules/align'
}

if (params.file == 'ubam') {
    include { bam2fastq  } from '../modules/align'
}

include {align} from '../modules/align'

// WORKFLOW

workflow alignment {
    take:
    sample_channel

    main: 
    if (params.style == 'ont' && params.file == 'fastq') {
        fastq_ch = cat(sample_channel)
    }

    else if (params.style == 'ont' && params.file == 'ubam') {
        fastq_ch= bam2fastq(sample_channel)

    }
    else {
        fastq_ch = sample_channel
    }

    aligned_bam_bai_channel = align(fastq_ch)

    emit:
    aligned_bam_bai_channel
}