

process bcf_snv {
    publishDir params.output, mode:'copy'

    beforeScript 'module load bioinfo-tools bcftools'
    
    cpus 8
    time '5h'

    input:
    path(bam)

    output:
    path "${bam.baseName}.snv.vcf", emit: snvfile

    script:
    """
    bcftools mpileup -X ${params.style} --threads 4 -f ${params.ref} ${bam} --annotate FORMAT/AD,FORMAT/DP | bcftools call --threads 4 -mv -P0.01 -Ob -o ${bam.baseName}.snv.vcf
    """
}