# LOMPE
LOng-read Multi-omics PipelinE

nextflow version 21.10.6

Workflow: 

minimap2, samtools, sniffles 1, BCFtools, CNVpytor, nanopolish, VEP, vcftools, custom database annotation and filtering of variants

Includes custom analysis of methylation calls from ONT fast5 files

QC: 

FastQC, picard

# RUN
nextflow run main.nf -config < > --fastq_folder < > --fast5_folder < optional > --sample_id < > 
--style < ont OR pb > 
--output <  > -with-trace 

requires;

pacbio (pb) or nanopore (ont) fastq files

For methylation calling:

Basecalled fast5 files (ont)

# Installation
git clone

modify config to your needs (executor, reference files, databases)

run!