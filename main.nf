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
    //include { meth_find } from './modules/run_methylation'
    
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
    
    //SNV calling 
    bcf_snv(align.out.bamfile)
    snv_filter(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(snv_filter.out)

    //phasing
    phase_it(align.out.bamfile, align.out.baifile, annotate_snvs.out)

    //add methylation info to bam
    fast5_folder= Channel.fromPath("${params.fast5_folder}")
    meth_index( fast5_folder, cat.out.fastq_file )
    meth_polish(cat.out.fastq_folder, phase_it.out.phased_bam, phase_it.out.phased_bai, meth_index.out.polish_index, meth_index.out.polish_index_fai, meth_index.out.polish_index_gzi, meth_index.out.polish_index_readdb )
    //split_on_chr(meth_polish.out.methylation_tsv)
    //meth_find(split_on_chr.out)
    //cat_meth(meth_find.out.collect())

    //extract phased methylation
    extract_meth(phase_it.out.phased_bam)

    //SV calling
    sniff(phase_it.out.phased_bam)
    pytor(phase_it.out.phased_bam, phase_it.out.phased_bai, bcf_snv.out.snvfile)
    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile)
    run_vep(combine.out)
    query(run_vep.out)
    filter(query.out)

    //QC
    picard(align.out.bamfile)
    fastqc(fastq_file)
    


}

// PB workflow
workflow pb {
    take: ${params.fastq_folder}

    main:
    fastq_folder = Channel.fromPath("${params.fastq_folder}/*gz").collect()
    cat(fastq_folder)
    align(cat.out.fastq_file)

    //SNV calling
    bcf_snv(align.out.bamfile)
    snv_filter(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(snv_filter.out)

    //phasing
    phase_it(align.out.bam, align.out.bai, annotate_snvs.out)

    //SV calling
    sniff(align.out.bamfile)
    pytor(align.out.bamfile, align.out.baifile, bcf_snv.out.snvfile)
    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile )
    run_vep(combine.out)

    query(run_vep.out)
    filter(query.out)

    //QC
//    picard(align.out.bamfile)
 //   fastqc(cat.out.fastq_file)

    

    
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