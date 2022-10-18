#!/usr/bin/env nextflow

/*
filter using DB and SVDB
*/


process query {
    publishDir params.output, mode: 'copy'

    input:
    path(combined)

    output:
    path "${combined.baseName}.query.vcf"; emit: queried

    script:
    if (${params.style} == 'pb') {
        DB = ${params.pb_DB}
    }
    else {
        DB = ${params.ont_DB}
    }
    """
    svdb --query --query_vcf ${combined} --db ${params.${DB}} --overlap -1 > ${combined.baseName}.query.vcf
    """
}

process filter {
    publishDir params.output, mode: 'copy'

    input:
    path(queried)

    output:
    path "${query.baseName}.filtered.vcf"; emit: filtered

    script:
    """
    python ./scripts/filter_rank.py ${queried} ${query.baseName}.filtered.vcf
    """
    
}