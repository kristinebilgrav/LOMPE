#!/usr/bin/env nextflow

/*
WORKFLOW for SNV calling and filtering
*/

// MODULES

include { deepvar } from '../modules/deepvariant'

if (params.bcftools) {
    include { bcf_snv ; filter_snvs } from '../modules/bcftools'
}
include {annotate_snvs} from '../modules/annotate'



// WORKFLOW

workflow snv_calling {
    take: aligned_bam_bai_channel

    main:
    deepvar(aligned_bam_bai_channel) 
    //bcf_snv(aligned_bam_bai_channel)
    filter_snvs(deepvar.out.deepvar_vcf) 

    //SNV annotate 
    annotate_snvs(filter_snvs.out.snv_filtered)
    bam_snv_channel = aligned_bam_bai_channel.join(annotate_snvs.out) 

    emit:
    bam_snv_channel

}