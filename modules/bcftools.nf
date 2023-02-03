

process bcf_snv {
    publishDir params.output, mode:'copy'

    beforeScript 'module load bioinfo-tools bcftools'
    errorStrategy 'ignore'
    
    cpus 8
    time '22h'

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
    bcftools mpileup -X ${LR} --threads 4 -f ${params.ref} ${bam} --annotate FORMAT/AD,FORMAT/DP | bcftools call --threads 4 -cv -P0.01 -p 1 -Ob |bcftools filter --exclude "QUAL<=0" |bcftools filter --exclude 'GT="ref"' > ${bam.baseName}.snv.vcf 

    """
}

process filter_snvs {
    publishDir params.output, mode:'copy'
    beforeScript 'module load bioinfo-tools bcftools'
    errorStrategy 'ignore'

    cpus 4
    time '12h'

    input:
    path(snvfile)

    output:
    path "${snvfile.baseName}.snv.filter.vcf", emit: snv_filter

    script:
    """
    bcftools filter --exclude "FORMAT/AD[0:0] < 3" $2.vcf | bcftools filter --exclude "FORMAT/AD[0:1] < 3" | bcftools filter --exclude "FORMAT/DP > 45" |bcftools filter --exclude "GT!='het'" > ${bam.baseName}.snv.filter.vcf
    """
}