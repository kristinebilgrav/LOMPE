#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process cat {
    tag "${params.style}:${SampleID}:cat" 

    input:

    tuple val(SampleID), file(fq_folder)
 

    output:
    tuple val(SampleID), file("${SampleID}.fastq.gz") , emit: fastq_file

    script:

    """ 
    zcat ${fq_folder}/fastq*/*gz > ${SampleID}.fastq
    gzip ${SampleID}.fastq
    """
}

process align {
    tag "${params.style}:${SampleID}:align"

    input:
    tuple val(SampleID), file(fastq)

    output:
    tuple val(SampleID), file( "${SampleID}.bam"), file( "${SampleID}.bam.bai")

    script: 
    """
    minimap2 -R '@RG\\tID:foo\\tSM:bar' -a -t ${task.cpus} --MD -x map-${params.style} ${params.ref} ${fastq} | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${SampleID}.bam && 
    samtools index -@ ${task.cpus} ${SampleID}.bam
    """
}