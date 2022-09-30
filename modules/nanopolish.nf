#!/usr/bin/env nextflow

/*
call methylation analysis using nanopolish
*/


process meth_polish {
    publishDir params.output, mode: 'copy'
    beforeScript 'module load hdf5'

    cpus 16
    time '42h'

    input:
    path(bam)
    path(fastq_gz)

    output:
    path , emit: methylation_tsv

    script:
    """
    export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin
    nanopolish call-methylation -r ${fastq_gz} -b ${bam} -g ${params.ref} --threads 16 > ${output}
    """
}

process split_on_chr {
    publishDir params.output, mode: 'copy'

    input: 
    path(methylation_tsv)

    output:
    path "chr*_meth.tsv", emit: methyl_chrs

    script:
    """
    for chr in $(seq 1 22) X Y
    do 
        head -n1 ${methylation_tsv} >chr$chr.txt
    done
    awk 'NR != 1 { print $0 >>("chr"$1".txt") }' ${methylation_tsv}
    """
}
