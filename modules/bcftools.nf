

process bcf_snv {
    tag "${params.style}:${SampleID}:bcftools"

    errorStrategy 'ignore'

    input:
    tuple val(SampleID), file(bam), file(bai)


    output:
    tuple val(SampleID), file("${bam.baseName}.snv.vcf"), emit: snvfile


    script:
    if (params.style == 'pb') 
        LR = 'pacbio-ccs'
    
    else 
        LR = 'ont'
    
    """
    bcftools mpileup -X ${LR} --threads 4 -f ${params.ref} ${bam} --annotate FORMAT/AD,FORMAT/DP | bcftools call --threads ${task.cpus} -cv -P0.01 -p 1 -Ob |bcftools filter --exclude "QUAL<=0" |bcftools filter --exclude 'GT="ref"' > ${bam.baseName}.snv.vcf 

    """
}

process filter_snvs {
    tag "${params.style}:${SampleID}:filter_snv"
    publishDir "${params.output}/${SampleID}_out/", mode:'copy'
    errorStrategy 'ignore'

    input:
    tuple val(SampleID), file(snvfile)

    output:
    tuple val(SampleID), file("${snvfile.baseName}.filter.vcf"), emit: snv_filtered

    script:
    """
    bcftools filter --exclude "FORMAT/AD[0:0] < 3" ${snvfile} | bcftools filter --exclude "FORMAT/AD[0:1] < 3" | bcftools filter --exclude "FORMAT/DP > 45" |bcftools filter --exclude "GT!='het'" > ${snvfile.baseName}.filter.vcf
    """
}