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
  vep -i ${queried} -o ${queried.baseName}.VEP.vcf ${params.vep_args} --vcf --fork 4 --per_gene --format vcf --no_stats
  """

}

process annotate_snvs {
  publishDir params.output, mode: 'copy'
  beforeScript 'module load bioinfo-tools vep'

  cpus 4
  time '3h'

  input:
  path(snv_filter)

  output: 
  path "${snv_filter.baseName}.snv.VEP.vcf", emit: annotated_snv_vcf

  script:
  """
  vep -i ${snv_filter} -o stdout --vcf  --per_gene --format vcf --no_stats --fork 4 --af --af_1kg --force_overwrite --af_gnomadg | filter_vep --filter "AF > 0.1 or gnomADg_AF > 0.1"  --force_overwrite -o ${snv_filter.baseName}.VEP.vcf
  """

}



