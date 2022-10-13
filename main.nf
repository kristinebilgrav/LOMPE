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
}

include { sniff } from './modules/sniffles'
include { run_vep } from './modules/annotate'
include { bcf_snv } from './modules/bcftools'
include { pytor } from './modules/pytor'
include { picard } from './modules/picard'
include { fastqc } from './modules/fastqc'
include {combine } from './modules/combine'

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
    
    run_vep(sniff.out)
    fastqc(bam)

}

workflow pb {
    
}

//methylation wf:
workflow meth {

    fastq = Channel.fromPath("${params.folder}")
    bam = Channel.fromPath(fastq.out.bamfile)
    bai = Channel.fromPath(fastq.out.baifile)
    meth_polish(fastq, bam)
    meth_find(meth_polish.out)

    combine(run_vep.out, pytor.out, meth_find.out)

}

// main workflow

workflow combine {
    
    combine(run_vep.out, pytor.out)
}


workflow LOMPE {
    call()
    if (params.style == 'ont') {
        meth(call.out)
    }
    if (params.style == 'pb') {
        combine(call.out)
    }


}
//completion handler
workflow.onComplete {
    println "Lompe complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check ${params.logfile}, maybe lompe fell apart :(")
}