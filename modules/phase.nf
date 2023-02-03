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
    path(annotated_snv_vcf)


    output:
    path "${bam.baseName}.phased.vcf.gz", emit: phased_vcf
    path "${bam.baseName}.phased.vcf.gz.tbi", emit: phased_vcf_tbi
    path "${bam.baseName}.haplotagged.bam", emit: phased_bam
    path "${bam.baseName}.haplotagged.bam.bai", emit: phased_bai


    script:
    """
    whatshap phase --tag=PS --ignore-read-groups -o ${bam.baseName}.phased.vcf --reference ${params.ref} ${bam} 
    bgzip ${bam.baseName}.phased.vcf
    tabix ${bam.baseName}.phased.vcf.gz
    whatshap haplotag -o ${bam.baseName}.haplotagged.bam --reference ${params.ref} ${bam.baseName}.phased.vcf.gz  ${bam} --output-threads=8
    samtools index -@ 8 ${bam.baseName}.haplotagged.bam
    """
    
}