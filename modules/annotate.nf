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
  path(queried)

  output:
  path "${queried.baseName}.VEP.vcf", emit: annotated_vcf

  script:
  """
  vep -i ${queried} -o ${queried.baseName}.VEP.vcf ${params.vep_args} 
  """

}



