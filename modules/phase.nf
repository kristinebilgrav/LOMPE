#!/usr/bin/env nextflow

/*
phase genome
*/

process phase_it {
    publishDir params.output, mode: 'copy'
    cpus 8
    time '16h'

    input:
    path(bam)
    path(snv_vcf)


    output:
    path "${bam.baseName}.phased.bam", emit: phased_bam
    path "${bam.baseName}.phased.bam.bai", emit: phased_bai


    script:
    if (params.style == 'pb') 
        phase_param = '/opt/margin_dir/params/phase/allParams.haplotag.pb-hifi.json'
    
    else 
        phase_param = '/opt/margin_dir/params/phase/allParams.haplotag.ont-r94g507.json'
    """
    margin phase ${bam} ${params.ref} ${snv_vcf} ${phase_param} -t 8 -o ${bam.baseName}.phased.bam
    samtools index -@ 8  ${bam.baseName}.phased.bam
    """
    
}