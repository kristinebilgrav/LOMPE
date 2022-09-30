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

include {align} from './modules/align.nf'
include { sniff } from './modules/sniffles.nf'
include { run_vep } from './modules/annotate.nf'
include { pytor } from './modules/pytor.nf'
include { picard } from './modules/picard.nf'
include { fastqc } from './modules/fastqc.nf'

// main workflow

workflow {
    align(params.folder)
    bam = Channel.fromPath('results/*bam')
    sniff(bam)
    pytor(bam)
    picard(bam)
    run_vep(sniff.out)
    fastqc(bam)
}