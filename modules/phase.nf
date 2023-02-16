#!/usr/bin/env nextflow

/*
phase genome
*/

process phase_it {
    tag "${params.style}:${SampleID}:phase"
    publishDir params.output, mode: 'copy'

    input:
    tuple val(SampleID), file(bamfile), file(baifile), file(annotated_snv_vcf)


    output:
    tuple val(SampleID),  file("${bam.baseName}.phased.vcf.gz"), emit: phased_vcf
    tuple val(SampleID), file("${bam.baseName}.phased.vcf.gz.tbi"), emit: phased_vcf_tbi
    tuple val(SampleID), file("${bam.baseName}.haplotagged.bam"), emit: phased_bam


    script:
    """
    whatshap phase --tag=PS --ignore-read-groups -o ${bam.baseName}.phased.vcf --reference ${params.ref} ${annotated_snv_vcf} ${bam} 
    bgzip ${bam.baseName}.phased.vcf
    tabix ${bam.baseName}.phased.vcf.gz
    whatshap haplotag -o ${bam.baseName}.haplotagged.bam --reference ${params.ref} ${bam.baseName}.phased.vcf.gz  ${bam} --output-threads=8
    """
    
}

process bamindex {
    tag "${params.style}:${SamplID}:index"
    publishDir params.output, mode: 'copy'

    input:
    tuple val(SampleID), file(bamfile)

    output:
    tuple val(SampleID), file(bamfile), file("${bam.baseName}.bam.bai")

    script:
    """
    samtools index -@ ${task.cpus} ${bamfile}
    """
}