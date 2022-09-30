#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process align {
    publishDir params.output, mode:'copy'
    beforeScript 'module load minimap2'

    cpus 16
    time '4h'

    input:
    path()

    output:
    path "${params.sample_id}.bam", emit: bamfile
    path "${params.sample_id}.bam.bai", emit: baifile

    script:
    """
    minimap2 -R '@RG\tID:foo\tSM:bar' -a -t 16 --MD -x map-${params.style} ${params.ref} ${}/* | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${params.sample_id}.bam && \
    samtools index ${params.sample_id}.bam
    """
}