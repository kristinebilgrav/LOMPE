#!/usr/bin/env nextflow

/*
main pipeline script
*/

nextflow.enable.dsl = 2



if ( params.input.endsWith('csv') ) {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ row ->  "${row.SamplePath}"  }
        .set{sample_channel}

}
else if (params.input.endsWith('bam')) {
    Channel
        .fromPath(params.input)
        .set{sample_channel}
        //begin at annotion
}

else if (new File(params.input).exists()){
    Channel
        .fromPath(params.input)
        .set{sample_channel}
}
else {
    println('input error, takes csv samplesheet, bam file or folder with gzipped fastq files')
}

log.info """\
LOng-read Multiomic PipelinE
----------------------------
input : ${params.input}
sample_id : ${params.sample_id}
output :  ${params.output}
style : ${params.style}

"""


// include modules
if (!params.input.endsWith('bam'))
    include {cat ; align} from './modules/align'

//include ont methylation annotation
if (params.style == 'ont') 
    include { meth_index ; meth_polish ; call_meth} from './modules/nanopolish'
    include { bamindex as index_methbam } from './modules/phase'
    
else if ( params.input.endsWith('bam'))
    include { cpg_tools } from './modules/cpg_tools'


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
    align(cat.out.fastq_file)
    
    //SNV calling 
    bcf_snv(align.out.bamfile)
    filter_snvs(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    phase_it(align.out.bamfile, align.out.baifile, annotate_snvs.out)
    bamindex(phase_it.out.phased_bam)

    //add methylation info to bam
    fast5_folder = "${sample_channel}/fast5*"
    meth_index( fast5_folder, cat.out.fastq_file )
    meth_polish(cat.out.fastq_file, phase_it.out.phased_bam, bamindex.out.bai, meth_index.out.polish_index, meth_index.out.polish_index_fai, meth_index.out.polish_index_gzi, meth_index.out.polish_index_readdb )
    index_methbam(meth_polish.out.methylation_bam)

    //extract phased methylation
    call_meth(meth_polish.out.methylation_bam, index_methbam.out.bai)


    //SV calling
    sniff(meth_polish.out.methylation_bam, index_methbam.out.bai)
    pytor(meth_polish.out.methylation_bam, index_methbam.out.bai, bcf_snv.out.snvfile)
    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile)
    run_vep(combine.out)
    query(run_vep.out)
    filter_query(query.out)

    //QC
    picard(align.out.bamfile)
    fastqc(cat.out.fastq_file)
    


}

// PB workflow
workflow pb_fastq {
    take: sample_channel

    main:
    //fastq_folder = Channel.fromPath("${params.fastq_folder}/*gz").collect()
    cat(sample_channel)
    align(cat.out.fastq_file)

    //SNV calling
    bcf_snv(align.out.bamfile)
    filter_snvs(bcf_snv.out)

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    phase_it(align.out.bam, align.out.bai, annotate_snvs.out)
    bamindex(phase_it.out.phased_bam)

    //SV calling
    sniff(phase_it.out.phased_bam, bamindex.out)
    pytor(phase_it.out.phased_bam, bamindex.out.bai, bcf_snv.out.snvfile)
    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile )
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
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    bai = "${sample_channel}.bai"
    phase_it(sample_channel, bai , annotate_snvs.out)
    bamindex(phase_it.out.phased_bam)

    //extract methylation from bam
    cpg_tools(phase_it.out.phased_bam, bamindex.out.bai)


    //SV calling
    sniff(phase_it.out.phased_bam, bamindex.out.bai)
    pytor(phase_it.out.phased_bam, bamindex.out.bai, bcf_snv.out.snvfile)
    combine(sniff.out.sniff_vcf, pytor.out.pytor_vcffile )
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
        ont(sample_channel)
    //    ont("${params.fastq_folder}")
    
    if (params.style == 'pb') 
        if (params.input.endsWith('bam')){
            println('starting pb bam channel')
            pb_bam(sample_channel)
        }
        else {
            pb_fastq(sample_channel)
        }
    
}

//completion handler
workflow.onComplete {
    println "LOMPE complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check ${params.logfile}, maybe lompe fell apart :(")
}