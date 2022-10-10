#!/usr/bin/env nextflow

/*
main pipeline script
*/

nextflow.enable.dsl = 2


log.info """\
Long-read pipeline
------------------
config: ${params.config}
folder(s) : ${params.folder}
sample_id : ${params.sample_id}
output :  ${params.output}

"""



// include modules

include {align} from './modules/align'
include { meth_polish } from './modules/nanopolish'
include { sniff } from './modules/sniffles'
include { run_vep } from './modules/annotate'
include {bcf_snv } from './modules/bcftools'
include { pytor } from './modules/pytor'
include { picard } from './modules/picard'
include { fastqc } from './modules/fastqc'
include {combine } from './modules/combine'

// main workflow

workflow {
    fastq = Channel.fromPath("${params.folder}")
    align(fastq)
    bam = Channel.fromPath('results/*bam')
    meth_polish(fastq, bam)
    meth_find(meth_polish.out)
    sniff(bam)
    bcf_snv(bam)
    pytor(bam, bcf_snv.out.snvfile)
    
    picard(bam)
    
    run_vep(sniff.out)
    fastqc(bam)

    combine(run_vep.out, pytor.out, meth_find.out)
}