#!/usr/bin/env nextflow

/*
identify SNVs with deepvariant
*/

process deepvar {
    tag "${params.style}:${SampleID}:deepvariant"
    publishDir "${params.output}/${SampleID}_out/", mode: 'copy'

    input:
    tuple val(SampleID), file(bamfile)
    

    output:
    tuple val(SampleID), file("${SampleID}_deepvariant.vcf"), emit: deepvar_vcf
    tuple val(SampleID), file("${SampleID}_deepvariant_g.vcf.gz"), emit: deepvar_gvcf


    script:
        if (params.style == 'pb') 
        LR = 'PACBIO'
    
    else 
        LR = 'ONT_R104'
    """
    /opt/deepvariant/bin/run_deepvariant --model_type=${LR} --ref=${params.ref} --reads=${bamfile} --output_vcf=${SampleID}_deepvariant.vcf.gz --output_gvcf=${SampleID}_deepvariant_g.vcf.gz --num_shards=${tasks.cpu}
    """
}
