#!/usr/bin/env nextflow

/*
extract methylation info from ONT using modbam2bed
*/

process modbam2bed {
    tag "${params.style}:${SampleID}:methbam2bed"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(bam), file(bai), file(fai)

    output:
    tuple val(SampleID), file("${SampleID}.methpos.bed"), emit: methpos_bed

    script:
    """
    python ${params.LOMPE_home}/scripts/modbam2bed.py ${bam} ${fai} > ${SampleID}.methpos.bed
    """

}