#!/usr/bin/env nextflow

/*
combine variant calls
*/

process combine {
    input:
    path(annotated_vcf)
    path(pytor)
    path(methyl)

    output:
    path "${annotated_vcf.simpleName}.output.vcf", emit: combined

    script:
    """
    cat 
    sort
    vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz
    """

}