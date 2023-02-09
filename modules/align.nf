#!/usr/bin/env nextflow

/*
alignment of fastq files using minimap2
*/

process cat {


    cpus 2
    time '2h'

    input:
    //path fq_files
    path(fq_folder)

    output:
    path "${sample_id}.fastq.gz", emit: fastq_file

    script:
    String path = fq_folder 
    sample_id = path.tokenize('/')[-1]
    """ 
    zcat ${fq_folder}/fastq_pass/*gz > ${sample_id}.fastq
    gzip ${sample_id}.fastq
    """
}

process align {


    input:
    path(fastq)

    output:
    path "${fastq.baseName}.bam", emit: bamfile
    path "${fastq.baseName}.bam.bai", emit: baifile

    script: 
    """
    minimap2 -R '@RG\\tID:foo\\tSM:bar' -a -t ${task.cpus} --MD -x map-${params.style} ${params.ref} ${fastq} | samtools view -Sbh - | samtools sort -m 4G -@16 - > ${fastq.baseName}.bam && 
    samtools index -@ ${task.cpus} ${fastq.baseName}.bam
    """
}