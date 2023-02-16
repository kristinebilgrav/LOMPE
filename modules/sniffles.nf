#!/usr/bin/env nextflow

/*
call SVs using sniffles 
*/

process sniff {
    tag "${params.style}:${SampleID}:sniffles"

    input:
    tuple val(SampleID), file(bam), file(bai)

    output:
    tuple val(SampleID), file("${bam.baseName}.sniffles.vcf"), emit: sniff_vcf


    shell:
    """
    sniffles -m ${bam} -v ${bam.baseName}.sniffles.vcf -l 100 -t ${task.cpus} -s 3 --genotype --cluster 

    """
}