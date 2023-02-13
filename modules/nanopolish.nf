#!/usr/bin/env nextflow

/*
call methylation analysis using nanopolish
*/

process meth_index {
    tag "${params.style}:${params.sample_id}:nanopolish_index"
    publishDir params.output, mode: 'copy'

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
    nanopolish index -d ${fast5_folder}/fast5_pass/ ${fastq_file}
    """
}


process meth_polish {
    tag "${params.style}:${params.sample_id}:nanopolish_call(bam)"
    publishDir params.output, mode: 'copy'

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
    nanopolish call-methylation -r ${fastq_file} -b ${bam} -g ${params.ref} --threads ${task.cpus} --methylation cpg --modbam-output-name ${bam.baseName}.methyl.bam
    """
}

process call_meth {
    tag "${params.style}:${params.sample_id}:nanopolish_call(tsv)"
    publishDir params.output, mode: 'copy'

    input:
    path(fastq_file)
    path(bam)
    path(bai)
    path(polish_index)
    path(polish_index_fai)
    path(polish_index_gzi)
    path(polish_index_readdb)

    output:
    path "${bam.baseName}.methylsites.tsv", emit: methylation_tsv
    path "${bam.baseName}.methylfrequency.tsv", emit: methylation_frequency

    script:
    """
    export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin
    nanopolish call-methylation -r ${fastq_file} -b ${bam} -g ${params.ref} --threads ${task.cpus} >  ${bam.baseName}.methylsites.tsv
    scripts/calculate_methylation_frequency.py ${bam.baseName}.methylsites.tsv > ${bam.baseName}.methylfrequency.tsv
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
