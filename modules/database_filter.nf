#!/usr/bin/env nextflow

/*
filter using DB and SVDB
*/


process query {
    tag "${params.style}:${SampleID}:SVDB"

    input:
    tuple val(SampleID), file(combined)

    output:
    tuple val(SampleID), file("${combined.baseName}.query.vcf")

    script:
    if (params.style == 'pb') 
        DB = params.pb_DB

    else 
        DB = params.ont_DB

    """
    svdb --query --query_vcf ${combined} --db ${DB} --overlap -1 > ${combined.baseName}.query.vcf
    """
}

process filter_query {
    tag "${params.style}:${SampleID}:filter"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(queried)

    output:
    tuple val(SampleID), file("${queried.baseName}.filtered.vcf"), emit: filtered

    script:
    """
    python ${params.LOMPE_home}/scripts/filter_rank.py ${queried} ${queried.baseName}.filtered.vcf
    """
    
}