#!/usr/bin/env nextflow

/*
combine variant calls
*/

process sort_zip {
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

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
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(sniff_vcf), file(pytor)



    output:
    tuple val(SampleID),  file("${sniff_vcf.simpleName}.output.vcf")

    script:
    """
    svdb --merge --vcf ${sniff_vcf} ${pytor} > ${sniff_vcf.simpleName}.output.vcf
    """

}