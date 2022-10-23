

process bcf_snv {
    publishDir params.output, mode:'copy'

    beforeScript 'module load bioinfo-tools bcftools'
    errorStrategy 'ignore'
    
    cpus 8
    time '10h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.snv.vcf", emit: snvfile

    script:
    if (params.style == 'pb') 
        LR = 'pacbio-ccs'
    
    else 
        LR = 'ont'
    
    """
    bcftools mpileup -X ${LR} --threads 4 -f ${params.ref} ${bam} --annotate FORMAT/AD,FORMAT/DP | bcftools call --threads 4 -mv -P0.01 -Ob -o ${bam.baseName}.snv.vcf
    """
}