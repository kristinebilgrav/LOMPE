#!/usr/bin/env nextflow

/*
call CNVs using pytor
*/

process pytor {
    publishDir params.output, mode:'copy'

    beforeScript 'module load bioinfo-tools bcftools'
    
    cpus 8
    time '5h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.snv.vcf", emit: snvfile
    path "${bam.baseName}.pytor", emit: pytorfile

    script:
    """

    bcftools mpileup -X ${params.style} --threads 4 -f ${params.ref} ${bam} --annotate FORMAT/AD,FORMAT/DP | bcftools call --threads 4 -mv -P0.01 -Ob -o ${bam}.snv.vcf

    cnvpytor -root ${bam.baseName}.pytor -rd ${bam.baseName}
    cnvpytor -root ${bam.baseName}.pytor -gc ${params.ref}
    cnvpytor -root ${bam.baseName}.pytor -his 20000 200000
    cnvpytor -root ${bam.baseName}.pytor -partition 20000 200000
    cnvpytor -root ${bam.baseName}.pytor -call 20000 > ${bam.baseName}.pytor
    cnvpytor -root ${bam.baseName}.pytor -snp ${bam.baseName}.snv.vcf -nofilter
    cnvpytor -root ${bam.baseName}.pytor -baf 20000 200000
    cnvpytor -root ${bam.baseName}.pytor -call combined 20000 > ${bam.baseName}.pytor.combined.out

    """
}