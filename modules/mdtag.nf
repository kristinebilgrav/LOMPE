#!/usr/bin/env nextflow

/*
MD tag genome
*/



process mdtag {
    tag "${params.style}:${SampleID}:MDtag"

    input:
    tuple val(SampleID), file(bam), file(bai)


    output:
    tuple val(SampleID), file("${bam.baseName}.md.bam"), file("${bam.baseName}.md.bam.bai")


    script:
    """
    samtools calmd -@ ${task.cpus} -b ${bam} --reference ${params.ref} > ${bam.baseName}.md.bam
    samtools index -@ ${task.cpus} ${bam.baseName}.md.bam
    """

}
