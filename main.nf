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

if (params.style == 'ont') {
    include { index ; meth_polish } from './modules/nanopolish'
    include { meth_find } from './modules/run_methylation'
    include { combine_ont } from './modules/combine'
}
else {
    include { combine_pb } from './modules/combine'
}

include { sniff } from './modules/sniffles'
include { run_vep } from './modules/annotate'
include { bcf_snv } from './modules/bcftools'
include { pytor } from './modules/pytor'
include { picard } from './modules/picard'
include { fastqc } from './modules/fastqc'


//workflows
workflow call {
    fastq = Channel.fromPath("${params.folder}/*fastq.gz")
    align(fastq)

    bam = Channel.fromPath(fastq.out.bamfile)
    bai = Channel.fromPath(fastq.out.baifile)

    sniff(bam)
    bcf_snv(bam)
    pytor(bam, bcf_snv.out.snvfile)

    picard(bam)
    fastqc(bam)

}



//methylation wf:
workflow ont {

    fastq = Channel.fromPath("${params.fastq_folder}")
    fast5= Channel.fromPath("${params.fast5_folder}")
    meth_polish(fastq, fast5)

    bam = Channel.fromPath(fastq.out.bamfile)
    bai = Channel.fromPath(fastq.out.baifile)
    meth_find(meth_polish.out, bam, bai )

    combine_ont(sniff.out, pytor.out, meth_find.out)
    run_vep(combine_ont.out)

}

// main workflow

workflow pb {
    
    combine_pb(sniff.out, pytor.out)
    run_vep(combine_pb.out)
}


workflow LOMPE {
    call(params.fastq_folder)
    if (params.style == 'ont') {
        ont(call.out)
    }
    if (params.style == 'pb') {
        pb(call.out)
    }

}

//completion handler
workflow.onComplete {
    println "LOMPE complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check ${params.logfile}, maybe lompe fell apart :(")
}