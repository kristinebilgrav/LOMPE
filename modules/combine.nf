#!/usr/bin/env nextflow

/*
combine variant calls
*/

process sort_zip {
    publishDir params.output, mode: 'copy'

    input:
    path(vcf)

    output:
    path "${vcf.baseName}.sort.vcf.gz", emit: vcffile
    path  "${vcf.baseName}.sort.vcf.gz.tbi", emit: vcffile_tbi

    script:
    """
    vcf-sort ${vcf} 2> /dev/null > ${vcf.baseName}.sort.vcf
    bgzip ${vcf.baseName}.sort.vcf
    tabix -p vcf ${vcf.baseName}.sort.vcf.gz
    """

}



process combine {
    tag "${params.style}:${SampleID}:SVDB"

    input:
    tuple val(SampleID), file(sniff_vcf), file(pytor)



    output:
    tuple val(SampleID),  file("${sniff_vcf.simpleName}.output.vcf"), emit: combined

    script:
    """
    svdb --merge --vcf ${sniff_vcf} ${pytor} > ${sniff_vcf.simpleName}.output.vcf
    """

}