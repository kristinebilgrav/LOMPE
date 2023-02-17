#!/usr/bin/env nextflow

/*
output methylation from bamfile to readable file
*/

process cpg_tools  {
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'
    tag "${params.style}:${SampleID}:cpg-tools"

    input:
    tuple val(SampleID), file(bam), file(bai)


    output:
    tuple val(SampleID),  file("${bam.baseName}*bed"), file("${bam.baseName}.*bw"), emit: methyl_sites

    script:
    """
    /miniconda/bin/python /pb-CpG-tools/aligned_bam_to_cpg_scores.py  -b ${bam} -f ${params.ref} -o ${bam.baseName} -t ${task.cpus} -q 10 -p model -d /pb-CpG-tools/pileup_calling_model
    """

}