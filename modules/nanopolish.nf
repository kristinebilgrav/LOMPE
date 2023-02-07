#!/usr/bin/env nextflow

/*
call methylation analysis using nanopolish
*/

process meth_index {
    publishDir params.output, mode: 'copy'
    cpus 6
    time '42h'
    container = 'quay.io/biocontainers/nanopolish'

    input:
    path(fast5_folder)
    path(fastq_file)
    

    output:
    path "${fastq_file}.index", emit: polish_index
    path "${fastq_file}.index.fai", emit: polish_index_fai
    path "${fastq_file}.index.gzi", emit: polish_index_gzi
    path "${fastq_file}.index.readdb", emit: polish_index_readdb
    

    script:
    """
    nanopolish index -d ${fast5_folder}/ ${fastq_file}
    """
}


process meth_polish {
    publishDir params.output, mode: 'copy'
    beforeScript 'module load hdf5'
    container = 'quay.io/biocontainers/nanopolish'

    cpus 16
    time '42h'

    input:
    path(fastq_file)
    path(bam)
    path(bai)
    path(polish_index)
    path(polish_index_fai)
    path(polish_index_gzi)
    path(polish_index_readdb)

    output:
    path "${bam.baseName}.methyl.bam", emit: methylation_bam

    script:
    """
    export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin
    nanopolish call-methylation -r ${fastq_file} -b ${bam} -g ${params.ref} --threads 16 --methylation cpg --modbam-output-name ${bam.baseName}.methyl.bam
    """
}

process split_on_chr {

    input: 
    path(methylation_tsv)

    output:
    path "chr*_meth.tsv", emit: methyl_chrs

    script:
    """
    for chr in \$(seq 1 22) X Y
    do 
        head -n1 ${methylation_tsv} > chr\$chr.txt
    done
    awk 'NR != 1 { print \$0 >>("chr"\$1".txt") }' ${methylation_tsv}
    """
}
