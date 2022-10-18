#!/usr/bin/env nextflow

/*
combine variant calls
*/

process combine_ont {
    publishDir params.output, mode: 'copy'
    
    input:
    path(sniff)
    path(pytor)
    path(methyl)

    output:
    path "${sniff.simpleName}.output.vcf", emit: combined

    script:
    """
    vcf-merge ${sniff} ${pytor} ${methyl} | bgzip -c > ${annotated_vcf.simpleName}.output.vcf
    """

}

process combine_pb {
    publishDir params.output, mode: 'copy'

    input:
    path(annotated_vcf)
    path(pytor)

    output:
    path "${annotated_vcf.simpleName}.output.vcf", emit: combined

    script:
    """
    vcf-merge ${annotated_vcf} ${pytor}| bgzip -c > ${annotated_vcf.simpleName}.output.vcf
    """

}