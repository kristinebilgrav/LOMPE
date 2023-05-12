#!/usr/bin/env nextflow

/*
call methylation analysis using nanopolish
*/

process meth_index {
    tag "${params.style}:${SampleID}:nanopolish_index"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), val(samplefolder), file(fastq_file)
    

    output:
    tuple val(SampleID), file("${fastq_file}.index"), file("${fastq_file}.index.fai"), file("${fastq_file}.index.gzi"), file("${fastq_file}.index.readdb")


    script:
    """
    nanopolish index -d ${samplefolder}/fast5_pass/ ${fastq_file}
    """
}


process meth_polish {
    tag "${params.style}:${SampleID}:nanopolish_call(bam)"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(fastq_file), file(bam), file(bai), file(polish_index), file(polish_index_fai), file(polish_index_gzi), file(polish_index_readdb)

    output:
    tuple val(SampleID),  file("${bam.baseName}.methyl.bam"), emit: methylation_bam

    script:
    """
    nanopolish call-methylation -r ${fastq_file} -b ${bam} -g ${params.ref} --threads ${task.cpus} --methylation cpg --modbam-output-name ${bam.baseName}.methyl.bam
    """
}

process call_meth {
    tag "${params.style}:${SampleID}:nanopolish_call(tsv)"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(fastq_file), file(bam), file(bai), file(polish_index), file(polish_index_fai), file(polish_index_gzi), file(polish_index_readdb)

    output:
    tuple val(SampleID), file("${bam.baseName}.methylsites.tsv"), file("${bam.baseName}.methylfrequency.tsv")

    script:
    """
    export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin
    nanopolish call-methylation -r ${fastq_file} -b ${bam} -g ${params.ref} --threads ${task.cpus} >  ${bam.baseName}.methylsites.tsv
    nanopolish/scripts/calculate_methylation_frequency.py ${bam.baseName}.methylsites.tsv > ${bam.baseName}.methylfrequency.tsv
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
