#!/usr/bin/env nextflow

/*
analyze methylation using script
*/

process meth_find {
    publishDir params.output, mode: 'copy'
    cpus 1
    time '2h'

    input:
    path(methyl_chrs)

    output:
    path "${methyl_chrs.baseName}.methyl_pos.tsv", emit: methylated_chr

    script:
    """
    python ./script/methylation.py ${methyl_chrs} ${methyl_chrs.baseName}.methyl_pos.tsv
    """
}

process cat_meth {
    publishDir params.output, mode: 'copy'
    cpus 1
    time '1h'

    input:
    path(methylated_chr) //all file in one?

    output:
    path "${params.sample_id}.methylated.tsv", emit: methylation_allchr

    script:
    """
    cat ${methylated_chr} > ${params.sample_id}.methylated.tsv
    """
}