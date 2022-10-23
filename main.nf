#!/usr/bin/env nextflow

/*
main pipeline script
*/

nextflow.enable.dsl = 2


log.info """\
LOng-read Multiomic PipelinE
----------------------------
config: ${params.config}
folder(s) : ${params.fastq_folder}
sample_id : ${params.sample_id}
output :  ${params.output}
style : ${params.style}

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
include { query ; filter} from './modules/database_filter'

//workflows

//ONT wf, with methylation:
workflow ont {
    take: ${params.fastq_folder}

    main:
    fastq_file = Channel.fromPath("${params.fastq_folder}/*fastq.gz")
    align(fastq_file)

    sniff(align.out.bamfile)
    bcf_snv(align.out.bamfile)
    pytor(align.out.bamfile, align.out.baifile, bcf_snv.out.snvfile)

    picard(align.out.bamfile)
    fastqc(fastq_file)

    fastq_folder = Channel.fromPath("${params.fastq_folder}")
    fast5_folder= Channel.fromPath("${params.fast5_folder}")
    meth_polish(fastq_folder, fast5_folder)

    meth_find(meth_polish.out)

    combine_ont(sniff.out, pytor.out, meth_find.out)
    run_vep(combine_ont.out)

    query(run_vep.out)
    filter(query.out)


}

// PB workflow
workflow pb {
    take: ${params.fastq_folder}

    main:
    fastq_file = Channel.fromPath("${params.fastq_folder}/*fastq.gz")
    align(fastq_file)


    sniff(align.out.bamfile)
    bcf_snv(align.out.bamfile)
    pytor(align.out.bamfile, align.out.baifile, bcf_snv.out.snvfile)

    picard(align.out.bamfile)
    fastqc(fastq_file)

    combine_pb(sniff.out, pytor.out.pytor_vcffile)
    run_vep(combine_pb.out)

    query(run_vep.out)
    filter(query.out)

    emit:
    fastqc.out.QC


}

//main workflow
workflow {

    main:
    if (params.style == 'ont') 
        ont("${params.fastq_folder}")
    
    if (params.style == 'pb') 
        pb("${params.fastq_folder}")
   
    emit:
    filtered

}

//completion handler
workflow.onComplete {
    println "LOMPE complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check ${params.logfile}, maybe lompe fell apart :(")
}