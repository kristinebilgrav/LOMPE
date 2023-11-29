#!/usr/bin/env nextflow

/*
phase genome
*/

process phase_it {
    tag "${params.style}:${SampleID}:phase"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(bam), file(bai), file(annotated_snv_vcf)


    output:
    tuple val(SampleID),  file("${bam.baseName}_snvs_phased.vcf.gz"), emit: phased_vcf
    tuple val(SampleID), file("${bam.baseName}_snvs_phased.vcf.gz.tbi"), emit: phased_vcf_tbi
    tuple val(SampleID), file("${bam.baseName}.haplotagged.bam"), emit: phased_bam


    script:
    """
    whatshap phase --tag=PS -o ${bam.baseName}_snvs_phased.vcf --reference ${params.ref} ${annotated_snv_vcf} ${bam} 
    bgzip ${bam.baseName}_snvs_phased.vcf
    tabix ${bam.baseName}_snvs_phased.vcf.gz
    whatshap haplotag -o ${bam.baseName}.haplotagged.bam --reference ${params.ref} ${bam.baseName}_snvs_phased.vcf.gz  ${bam} --output-threads=${task.cpus}
    """
    
}

process bamindex {
    tag "${params.style}:${SampleID}:index"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(bamfile)

    output:
    tuple val(SampleID), file(bamfile), file("${bamfile}.bai")

    script:
    """
    samtools index -@ ${task.cpus} ${bamfile}
    """
}