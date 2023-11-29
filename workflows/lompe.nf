#!/usr/bin/env nextflow

/*
LOMPE 
*/

nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



//if (params.style != 'pb' ||  'ont' ) { exit 1, 'style must be either "pb" or "ont" '}
//if (params.file != 'bam' ||  'fastq' || 'ubam') { exit 1, 'file must be either bam, ubam or fastq'}
if ( params.input.endsWith('csv') ) { 
    /*
    if samplesheet - need to end with csv. 
    treat differently if its pacbio (pb) bam file (do not need alignment)
    otherwise SampleID and SamplePath saved to sample_channel or 
    */
    if (params.file == 'bam' ||  'ubam' ) {
        Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map{ row ->  tuple(row.SampleID, file(row.SamplePath), file(row.SamplePath.tokenize('.')[0]+'*.bai'))  }
            .set{sample_channel}
    }
    else {
        Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ row ->  tuple(row.SampleID, file(row.SamplePath))  }
        .set{sample_channel}
    }
    
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def mandatoryParams = [
    "input",
    "style",
    "ref",
    "file",
    "output"
]


def missingParamsCount = 0
for (p in mandatoryParams.unique()) {
    if (params[p] == null) {
        println("params." + p + " not set.")
        missingParamsCount += 1
    }
}

if (missingParamsCount>0) {
    error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES




//include ont methylation annotation
if (params.file == 'fastq' && params.style == 'ont') {
    include { meth_index ; meth_polish ; call_meth} from '../modules/nanopolish'
    include { bamindex as index_methbam } from '../modules/phase'
}    
else if ( params.style == 'pb'){
    include { cpg_tools } from '../modules/cpg_tools'
    include { mdtag } from '../modules/mdtag'
}

include { combine } from '../modules/combine'



include { pytor } from '../modules/pytor'

include { query ; filter_query} from '../modules/database_filter'


// SUBWORKFLOWS

include { alignment } from '../subworkflows/wf-alignment'
include { snv_calling } from '../subworkflows/wf-SNV'
include { phasing } from '../subworkflows/wf-phase'
include { sv_calling } from '../subworkflows/wf-SVs'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOMPE {

    if (params.file == 'bam'){
        println('Starting from BAM, only variant calling')
        main:
        snv_calling(sample_channel)
        sv_calling(sample_channel)
    }

    else {
        println('Starting with raw data')
        main:
        alignment(sample_channel)
        snv_calling(alignment.out)
        phasing(snv_calling.out)
        sv_calling(phasing.out)
        //methylation_calling()
    }       

}
    
    

   


log.info """\
LOng-read Multiomic PipelinE
----------------------------
input : ${params.input}
output :  ${params.output}
style : ${params.style}

"""

//completion handler
workflow.onComplete {
    println "LOMPE complete at $workflow.complete"
    log.info (workflow.success ? "Done! Lompe is filled with aligned and analyzed files at ${params.output}!" : "Fail. Check the trace file, maybe lompe fell apart :(")
}




//else if (params.style == 'pb') { }

//if (params.style == 'ont') {}