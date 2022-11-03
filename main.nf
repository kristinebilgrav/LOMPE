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

include {cat ; align} from './modules/align'

if (params.style == 'ont') {
    include { meth_index ; meth_polish } from './modules/nanopolish'
    include { meth_find } from './modules/run_methylation'
    
}

include { combine } from './modules/combine'
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
    fastq_folder = Channel.fromPath("${params.fastq_folder}/*").collect()
    cat(fastq_folder)
    align(cat.out.fastq_file)

    fast5_folder= Channel.fromPath("${params.fast5_folder}")
    meth_index( fast5_folder, cat.out.fastq_file )
    meth_polish(cat.out.fastq_folder, bam, bai, meth_index.out.polish_index, meth_index.out.polish_index_fai, meth_index.out.polish_index_gzi, meth_index.out.polish_index_readdb )
    meth_find(meth_polish.out.methylation_tsv)
    
    
    sniff(align.out.bamfile)
    bcf_snv(align.out.bamfile)
    pytor(align.out.bamfile, align.out.baifile, bcf_snv.out.snvfile)
    sort_zip(pytor.out)

    picard(align.out.bamfile)
    fastqc(fastq_file)

    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile)
    run_vep(combine_ont.out)

    query(run_vep.out)
    filter(query.out)


}

// PB workflow
workflow pb {
    take: ${params.fastq_folder}

    main:
    fastq_folder = Channel.fromPath("${params.fastq_folder}/*").collect()
    cat(fastq_folder)
    align(cat.out.fastq_file)


    sniff(align.out.bamfile)
    bcf_snv(align.out.bamfile)
    pytor(align.out.bamfile, align.out.baifile, bcf_snv.out.snvfile)
    sort_zip(pytor.out.pytor_vcffile)

    picard(align.out.bamfile)
    fastqc(fastq_file)

    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile )
    run_vep(combine_pb.out)

    query(run_vep.out)
    filter(query.out)

    
}

//main workflow
workflow {

    main:
    if (params.style == 'ont') 
        ont("${params.fastq_folder}")
    
    if (params.style == 'pb') 
        pb("${params.fastq_folder}")
   
    
}

//completion handler
workflow.onComplete {
    println "LOMPE complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check ${params.logfile}, maybe lompe fell apart :(")
}