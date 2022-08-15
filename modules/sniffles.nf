#!/usr/bin/env nextflow

/*
call SVs using sniffles 
*/

process sniff {
    publishDir params.outdir, mode:'copy'

    beforeScript 'module load bioinfo-tools Sniffles'

    input:
    path(bam)

    output:
    path "${bam.baseName}.sniffles.vcf"

    shell:
    """
    sniffles -m ${bam} -v ${bam.baseName}.sniffles.vcf -l 100 -t 16 -s 3 --genotype --cluster
    """
}