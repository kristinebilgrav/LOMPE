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
    path "${methyl_chrs.baseName}.fastQC.out/*"

    script:
    """
    python ./script/methylation.py chr_file.tsv ${output}
    """
}