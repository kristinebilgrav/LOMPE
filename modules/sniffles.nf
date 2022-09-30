#!/usr/bin/env nextflow

/*
call SVs using sniffles 
*/

process sniff {
    publishDir params.output, mode:'copy'

    beforeScript 'module load bioinfo-tools Sniffles'

    input:
    path(bam)

    output:
    path "${bam.baseName}.sniffles.vcf", emit: sniff_vcf

    shell:
    """
    sniffles -m ${bam} -v ${bam.baseName}.sniffles.vcf -l 100 -t 16 -s 3 --genotype --cluster
    """
}