#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process align {
    publishDir params.outdir, mode:'copy'

    cpus 16
    time '4h'

    input:
    path(fqfolder)

    output:
    path "${params.sample_id}.bam", emit: bamfile
    path "${params.sample_id}.bam.bai" 

    script:
    """
    if params.style
    minimap2 -R '@RG\tID:foo\tSM:bar' -a -t 16 --MD -x map-${params.style} ${params.ref} ${fqfolder}/* | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${params.sample_id}.bam |
    samtools index ${params.sample_id}.bam
    """
}