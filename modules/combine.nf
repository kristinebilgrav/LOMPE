#!/usr/bin/env nextflow

/*
combine variant calls
*/

process combine_ont {
    publishDir params.output, mode: 'copy'
    
    input:
    path(sniff_vcf)
    path(pytor)
    path(methyl)

    output:
    path "${sniff_vcf.simpleName}.output.vcf", emit: combined

    script:
    """
    vcf-merge ${sniff_vcf} ${pytor} ${methyl} | bgzip -c > ${sniff_vcf.simpleName}.output.vcf
    """

}

process combine_pb {
    publishDir params.output, mode: 'copy'

    input:
    path(sniff_vcf)
    path(pytor)

    output:
    path "${sniff_vcf.simpleName}.output.vcf", emit: combined

    script:
    """
    vcf-merge ${sniff_vcf} ${pytor}| bgzip -c > ${sniff_vcf.simpleName}.output.vcf
    """

}