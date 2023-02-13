#!/usr/bin/env nextflow

/*
phase genome
*/

process phase_it {
    tag "${params.style}:${params.sample_id}:phase"
    publishDir params.output, mode: 'copy'

    input:
    path(bam)
    path(baifile)
    path(annotated_snv_vcf)


    output:
    path "${bam.baseName}.phased.vcf.gz", emit: phased_vcf
    path "${bam.baseName}.phased.vcf.gz.tbi", emit: phased_vcf_tbi
    path "${bam.baseName}.haplotagged.bam", emit: phased_bam


    script:
    """
    whatshap phase --tag=PS --ignore-read-groups -o ${bam.baseName}.phased.vcf --reference ${params.ref} ${annotated_snv_vcf} ${bam} 
    bgzip ${bam.baseName}.phased.vcf
    tabix ${bam.baseName}.phased.vcf.gz
    whatshap haplotag -o ${bam.baseName}.haplotagged.bam --reference ${params.ref} ${bam.baseName}.phased.vcf.gz  ${bam} --output-threads=8
    """
    
}

process bamindex {
    tag "${params.style}:${params.sample_id}:index"
    publishDir params.output, mode: 'copy'
    cpus 8
    time '2h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.bam.bai", emit: bai

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}