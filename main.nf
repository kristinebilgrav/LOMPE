#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    kristinebilgrav/LOMPE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/kristinebilgrav/LOMPE  
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { LOMPE } from './workflows/lompe'


// WORKFLOW: Run main kristinebilgrav/LOMPE analysis pipeline
workflow  {
    LOMPE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//


/*
workflows:
two for PB and two for ONT, 
depending on fastq or bam file input


//ONT wf, fastq input with methylation calling from fast5:
workflow ont_fastq {
    take: sample_channel

    main:
    cat(sample_channel)
    aligned_bam_bai_channel = align(cat.out.fastq_file)
    
    //SNV calling 
    //bcf_snv(aligned_bam_bai_channel)
    deepvar(aligned_bam_bai_channel) 
    filter_snvs(deepvar.out.deepvar_vcf) 

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    bam_snv_channel = aligned_bam_bai_channel.join(annotate_snvs.out) //join
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out.phased_bam) //may have to join seperatly

    //add methylation info to bam
    index_channel = sample_channel.join(cat.out) //join
    index_files = meth_index( index_channel )
    fastq_bam = cat.out.join(phased_bam_bai_channel)
    call_methylation_channel = fastq_bam.join(index_files) //join 
    meth_polish(call_methylation_channel)
    methylation_bam_bai_channel = index_methbam(meth_polish.out.methylation_bam)

    //extract phased methylation
    call_meth(call_methylation_channel)
   
    //SV calling
    //join unless meth channel works
    sniff(methylation_bam_bai_channel)
    pytor_in_channel = methylation_bam_bai_channel.join(deepvar.out.deepvar_vcf)//join
    pytor(pytor_in_channel)
    combine_channel = sniff.out.join(pytor.out.pytor_vcf) //join
    combine(combine_channel)
    run_vep(combine.out)
    query(run_vep.out)
    filter_query(query.out)

    //QC
    picard(align.out)
    fastqc(cat.out)
    
}

//ONT wf, from uBAM:
workflow ont_bam {
    take: sample_channel

    main:
    bam2fastq(sample_channel)
    aligned_bam_bai_channel = align(bam2fastq.out.fastq_file)
    
    //SNV calling 
    //bcf_snv(aligned_bam_bai_channel)
    //filter_snvs(bcf_snv.out)
    deepvar(aligned_bam_bai_channel)
    filter_snvs(deepvar.out.deepvar_vcf)
    

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)

    //phasing
    bam_snv_channel = aligned_bam_bai_channel.join(annotate_snvs.out) //join
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out.phased_bam) //may have to join seperatly

    //extract phased methylation
    cpg_tools(phased_bam_bai_channel) //modkitbam??
   
    //SV calling
    //join unless meth channel works
    sniff(phased_bam_bai_channel)
    pytor_in_channel = phased_bam_bai_channel.join(deepvar.out.deepvar_vcf)//join
    pytor(pytor_in_channel)
    combine_channel = sniff.out.join(pytor.out.pytor_vcf) //join
    combine(combine_channel)
    run_vep(combine.out)
    query(run_vep.out)
    filter_query(query.out)

    //QC
    picard(align.out)
    fastqc(bam2fastq.out)
    
}


// PB workflow from fastq
workflow pb_fastq {
    take: sample_channel

    main:
    aligned_bam_bai_channel = align(sample_channel)
    
    //SNV calling 
    deepvar(aligned_bam_bai_channel)
    filter_snvs(deepvar.out.deepvar_vcf)
    //bcf_snv(mdtag.out)
    //filter_snvs(bcf_snv.out)

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
    combine_channel = sniff.out.join(pytor.out.pytor_vcf) //join
    combine(combine_channel )
    run_vep(combine.out)

    query(run_vep.out)
    filter_query(query.out)

    //QC
    picard(aligned_bam_bai_channel)
    fastqc(sample_channel)
    
}

//PB workflow from bam file with methylation marks
workflow pb_bam {
    take: sample_channel

    main:
    //SNV calling 
    mdtag(sample_channel)
    deepvar(aligned_bam_bai_channel)
    filter_snvs(deepvar.out.deepvar_vcf)
    //bcf_snv(mdtag.out)
    //filter_snvs(bcf_snv.out)
    

    //SNV annotate 
    annotate_snvs(filter_snvs.out)

    //phasing
    bam_snv_channel = mdtag.out.join(annotate_snvs.out) //join
    phase_it(bam_snv_channel)
    phased_bam_bai_channel = bamindex(phase_it.out.phased_bam) //may have to join seperatly

    //extract methylation from bam
    cpg_tools(phased_bam_bai_channel)

    //SV calling
    sniff(phased_bam_bai_channel)
    pytor_in_channel = phased_bam_bai_channel.join(deepvar.out.deepvar_vcf)//join
    pytor(pytor_in_channel)
    combine_channel = sniff.out.join(pytor.out.pytor_vcf) //join
    combine(combine_channel )
    run_vep(combine.out)

    query(run_vep.out)
    filter_query(query.out)

    //QC
    picard(phased_bam_bai_channel)

}


*/

