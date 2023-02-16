#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process cat {
    tag "${params.style}:${SampleID}:cat" 

    input:

    tuple val(SampleID), file(fq_folder)
 

    output:
    tuple val(SampleID), file("${params.sample_id}.fastq.gz") , emit: fastq_file

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
    tuple val(SampleID), file( "${fastq.baseName}.bam"), file( "${fastq.baseName}.bam.bai")

    script: 
    """
    minimap2 -R '@RG\\tID:foo\\tSM:bar' -a -t ${task.cpus} --MD -x map-${params.style} ${params.ref} ${fastq} | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${fastq.baseName}.bam && 
    samtools index -@ ${task.cpus} ${fastq.baseName}.bam
    """
}