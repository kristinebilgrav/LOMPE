#!/usr/bin/env nextflow

/*
main pipeline script
*/

nextflow.enable.dsl = 2



if ( params.input.endsWith('csv') ) { 
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ row ->  tuple(row.SampleID, file(row.SamplePath))  }
        .set{sample_channel}

}
else if (params.input.endsWith('bam')) {
    String path = params.input 
    SampleID = path.tokenize('/')[-1].tokenize('.')[0]
    Channel
        .fromPath(params.input)
        .set{tuple(SampleID, sample_channel)}
        //begin at annotion
}

else if (new File(params.input).exists()){

    String path = params.input 
    SampleID = path.tokenize('/')[-1]

    Channel
        .fromPath(params.input)
        .set{tuple(SampleID,sample_channel)}
}
else {
    println('input error, takes csv samplesheet, bam file or folder with gzipped fastq files')
}

//



// include modules
if (!params.input.endsWith('bam'))
    include {cat ; align} from './modules/align'

//include ont methylation annotation
if (params.style == 'ont') {
    include { meth_index ; meth_polish ; call_meth} from './modules/nanopolish'
    include { bamindex as index_methbam } from './modules/phase'
}    
else if ( params.style == 'pb'){
    include { cpg_tools } from './modules/cpg_tools'
}

include { combine } from './modules/combine'
include { sniff } from './modules/sniffles'
include { run_vep ; annotate_snvs} from './modules/annotate'
include { phase_it ; bamindex } from './modules/phase'
include { bcf_snv ; filter_snvs } from './modules/bcftools'
include { pytor } from './modules/pytor'
include { picard } from './modules/picard'
include { fastqc } from './modules/fastqc'
include { query ; filter_query} from './modules/database_filter'

//workflows

//ONT wf, with methylation:
workflow ont {
    take: sample_channel

    main:
    cat(sample_channel)
    aligned_bam_bai_channel = align(cat.out.fastq_file)
    
    //SNV calling 
    bcf_snv(aligned_bam_bai_channel)
    filter_snvs(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    bam_snv_channel = aligned_bam_bai_channel.join(annotate_snvs.out) //join
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out.phased_bam) //may have to join seperatly

    //add methylation info to bam
    index_channel = sample_channel.join(cat.out) //join
    index_files = meth_index( index_channel )
    call_methylation_channel = phased_bam_bai_channel.join(index_files) //join 
    meth_polish(call_methylation_channel)
    //meth_polish(cat.out.fastq_file, phase_it.out.phased_bam, bamindex.out.bai, meth_index.out.polish_index, meth_index.out.polish_index_fai, meth_index.out.polish_index_gzi, meth_index.out.polish_index_readdb )
    methylation_bam_bai_channel = index_methbam(meth_polish.out.methylation_bam)

    //extract phased methylation
    call_meth(call_methylation_channel)
    //call_meth(cat.out.fastq_file, meth_polish.out.methylation_bam, index_methbam.out.bai, meth_index.out.polish_index, meth_index.out.polish_index_fai, meth_index.out.polish_index_gzi, meth_index.out.polish_index_readdb)

    //SV calling
    //join unless meth channel works
    sniff(methylation_bam_bai_channel)
    pytor_in_channel = methylation_bam_bai_channel.join(bcf_snv.out)//join
    pytor(pytor_in_channel)
    combine_channel = sniff.out.join(pytor.out) //join
    combine(combine_channel)
    run_vep(combine.out)
    query(run_vep.out)
    filter_query(query.out)

    //QC
    picard(align.out)
    fastqc(cat.out.fastq_file)
    


}

// PB workflow
workflow pb_fastq {
    take: sample_channel

    main:
    cat(sample_channel)
    aligned_bam_bai_channel = align(cat.out.fastq_file)
    
    //SNV calling 
    bcf_snv(aligned_bam_bai_channel)
    filter_snvs(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    bam_snv_channel = aligned_bam_bai_channel.join(annotate_snvs.out) //join
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out.phased_bam) //may have to join seperatly

    //SV calling
    sniff(phased_bam_bai_channel)
    pytor_in_channel = phased_bam_bai_channel.join(bcf_snv.out)//join
    pytor(pytor_in_channel)
    combine_channel = sniff.out.join(pytor.out) //join
    combine(combine_channel )
    run_vep(combine.out)

    query(run_vep.out)
    filter_query(query.out)

    //QC
//    picard(align.out.bamfile)
 //   fastqc(cat.out.fastq_file)
    
}

workflow pb_bam {
    take: sample_channel

    main:
    //SNV calling 
    bcf_snv(sample_channel)
    filter_snvs(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(filter_snvs.out)

    //phasing
    bam_snv_channel = sample_channel.join(annotate_snvs.out) //join
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out) //may have to join seperatly

    //extract methylation from bam
    cpg_tools(phased_bam_bai_channel)

    //SV calling
    sniff(phased_bam_bai_channel)
    pytor_in_channel = phased_bam_bai_channel.join(bcf_snv.out)//join
    pytor(pytor_in_channel)
    combine_channel = sniff.out.join(pytor.out) //join
    combine(combine_channel )
    run_vep(combine.out)

    query(run_vep.out)
    filter_query(query.out)

    //QC
//    picard(sample_channel)

}

//main workflow
workflow {

    main:
    if (params.style == 'ont') 
        println('ONT workflow starting')
        sample_channel.view()
        ont(sample_channel)

    
    if (params.style == 'pb')
        
        sample_channel.view() 
        if (params.file == 'bam'){
            println('PB workflow starting from BAM')
            pb_bam(sample_channel)
        }
        else {
            println('PB workflow starting from fastq')
            pb_fastq(sample_channel)
        }
    
}


log.info """\
LOng-read Multiomic PipelinE
----------------------------
input : folder
sample_id : sample
output :  ${params.output}
style : ${params.style}

"""

//completion handler
workflow.onComplete {
    println "LOMPE complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check ${params.logfile}, maybe lompe fell apart :(")
}