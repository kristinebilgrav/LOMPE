#!/usr/bin/env nextflow

/*
call CNVs using pytor
*/

process pytor {
    tag "${params.style}:${SampleID}:CNVpytor"

    input:
    tuple val(SampleID), file(bam), file(bai), file(snvfile)

    output:
    tuple val(SampleID), file("${bam.baseName}.pytor.vcf"), file("${bam.baseName}.pytor.out") 

    script:
    """
    cnvpytor -root ${bam.baseName}.pytor -rd ${bam} 
    cnvpytor -root ${bam.baseName}.pytor -gc ${params.ref} 
    cnvpytor -root ${bam.baseName}.pytor -his 20000 200000 
    cnvpytor -root ${bam.baseName}.pytor -partition 20000 200000 
    cnvpytor -root ${bam.baseName}.pytor -call 20000 > ${bam.baseName}.pytor.out
    cnvpytor -root ${bam.baseName}.pytor -snp ${snvfile} -nofilter 
    cnvpytor -root ${bam.baseName}.pytor -baf 20000 200000 
    cnvpytor -root ${bam.baseName}.pytor -call combined 20000 > ${bam.baseName}.pytor.combined.out 

    
    cnvpytor -root ${bam.baseName}.pytor -view 20000 <<ENDL
    set print_filename ${bam.baseName}.pytor.vcf
    print calls
    ENDL

    """
}