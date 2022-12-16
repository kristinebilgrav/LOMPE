#!/usr/bin/env nextflow

/*
output methylation from bamfile to readable file
*/

process extract_meth {
    publishDir params.output, mode: 'copy'
    cpus 8
    time '16h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.bed", emit: methyl_sites
    path "${bam.baseName}.bw", emit: methyla_bwtrack

    script:
    """
    /miniconda/bin/python pb-CpG-tools/aligned_bam_to_cpg_scores.py  -b ${bam} -f ${params.ref} -o ${bam.baseName} -t 8 -d /pb-CpG-tools/pileup_calling_model
    """

}