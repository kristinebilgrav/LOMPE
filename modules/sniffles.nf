#!/usr/bin/env nextflow

/*
call SVs using sniffles 
*/

process sniff {
    beforeScript 'module load bioinfo-tools Sniffles'
    cpus 16
    time '5h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.sniffles.vcf", emit: sniff_vcf


    shell:
    """
    sniffles -m ${bam} -v ${bam.baseName}.sniffles.vcf -l 100 -t ${task.cpus} -s 3 --genotype --cluster 

    """
}