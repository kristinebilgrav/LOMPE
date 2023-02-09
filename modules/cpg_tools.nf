#!/usr/bin/env nextflow

/*
output methylation from bamfile to readable file
*/

process cpg_tools  {
    publishDir params.output, mode: 'copy'

    input:
    path(bam)
    path(bai)

    output:
    path "${bam.baseName}*bed", emit: methyl_sites
    path "${bam.baseName}.*bw", emit: methyl_bwtrack

    script:
    """
    /miniconda/bin/python /pb-CpG-tools/aligned_bam_to_cpg_scores.py  -b ${bam} -f ${params.ref} -o ${bam.baseName} -t ${task.cpus} -q 10 -p model -d /pb-CpG-tools/pileup_calling_model
    """

}