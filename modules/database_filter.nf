#!/usr/bin/env nextflow

/*
filter using DB and SVDB
*/


process query {
    tag "${params.style}:${params.sample_id}:SVDB"
    input:
    path(combined)

    output:
    path "${combined.baseName}.query.vcf"

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
    tag "${params.style}:${params.sample_id}:filter"
    publishDir params.output, mode: 'copy'

    input:
    path(queried)

    output:
    path "${queried.baseName}.filtered.vcf", emit: filtered

    script:
    """
    python ${params.LOMPE_home}/scripts/filter_rank.py ${queried} ${queried.baseName}.filtered.vcf
    """
    
}