#!/usr/bin/env nextflow

/*
annotate SV calls from sniffles
*/

process run_vep {
  publishDir params.output, mode: 'copy'
  beforeScript 'module load bioinfo-tools vep'

  cpus 4
  time '1h'

  input:
  path(sniff_vcf)

  output:
  path "${bam.baseName}.VEP.vcf", emit: annotated_vcf

  shell:
  """
  ${params.vep_path} -i ${} -o ${bam.baseName}.VEP.vcf ${params.vep_args} 
  """

}



