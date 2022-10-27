#!/usr/bin/env nextflow

/*
filter using DB and SVDB
*/


process query {
    publishDir params.output, mode: 'copy'

    input:
    path(combined)

    output:
    path "${combined.baseName}.query.vcf", emit: queried

    script:
    if (params.style == 'pb') 
        DB = params.pb_DB

    else 
        DB = params.ont_DB

    """
    svdb --query --query_vcf ${combined} --db ${DB} --overlap -1 > ${combined.baseName}.query.vcf
    """
}

process filter {
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