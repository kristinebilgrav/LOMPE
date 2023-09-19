#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process cat {
    tag "${params.style}:${SampleID}:cat" 
    publishDir "${params.output}/${SampleID}_out/", mode:'copy'

    input:

    tuple val(SampleID), file(fq_folder)
 

    output:
    tuple val(SampleID), file("${SampleID}.fastq.gz") , emit: fastq_file

    script:

    """ 
    zcat ${fq_folder}/*fastq* > ${SampleID}.fastq
    gzip ${SampleID}.fastq
    """
}

process bam2fastq {
    tag "${params.style}:${SampleID}:bam2fastq"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(bam_folder)


    output:
    tuple val(SampleID), file("${SampleID}.fastq.gz") , emit: fastq_file

    script:

    """
    samtools merge -o ${SampleID}.bam ${bam_folder}/*.bam
    samtools index -@ ${task.cpus} ${SampleID}.bam
    samtools fastq -TMM,ML ${SampleID}.bam > ${SampleID}.fastq
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
    //add Y flag for making supplementary alignments soft clips?
    """
    minimap2 -R '@RG\\tID:foo\\tSM:bar' -a -t ${task.cpus} --MD -x map-${params.style} -y ${params.ref} ${fastq} | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${SampleID}.bam && 
    samtools index -@ ${task.cpus} ${SampleID}.bam
    """
}