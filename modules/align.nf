#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process cat {


    cpus 2
    time '2h'

    input:
    path fq_files

    output:
    path "${params.sample_id}.fastq.gz", emit: fastq_file

    script:
    """
    zcat ${fq_files} > ${params.sample_id}.fastq  
    gzip ${params.sample_id}.fastq
    """
}

process align {


    input:
    path(fastq)

    output:
    path "${params.sample_id}.bam", emit: bamfile
    path "${params.sample_id}.bam.bai", emit: baifile

    script: 
    """
    minimap2 -R '@RG\\tID:foo\\tSM:bar' -a -t ${task.cpus} --MD -x map-${params.style} ${params.ref} ${fastq} | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${params.sample_id}.bam && 
    samtools index ${params.sample_id}.bam
    """
}