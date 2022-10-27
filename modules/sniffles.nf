#!/usr/bin/env nextflow

/*
call SVs using sniffles 
*/

process sniff {
    publishDir params.output, mode:'copy'
    beforeScript 'module load bioinfo-tools Sniffles'
    cpus 16
    time '5h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.sniffles.sort.vcf.gz", emit: sniff_vcf
    path "${bam.baseName}.sniffles.sort.vcf.gz.tbi", emit: sniff_vcf_tbi

    shell:
    """
    sniffles -m ${bam} -v ${bam.baseName}.sniffles.vcf -l 100 -t 16 -s 3 --genotype --cluster &&
    vcf-sort ${bam.baseName}.sniffles.vcf  2> /dev/null > ${bam.baseName}.sniffles.sort.vcf
    bgzip ${bam.baseName}.sniffles.sort.vcf
    tabix -p vcf ${bam.baseName}.sniffles.sort.vcf.gz

    """
}