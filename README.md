# LOMPE
LOng-read Multi-omics PipelinE

* Under development *

nextflow version 21.10.6 (version important)

Workflow: 
fastq -> bam (minimap2)
snv calling and annotation (bcftools)
phasing and haplotagging of bam file (whatshap)
methylation calling (ont) (nanopolish)
sv calling and cnv calling (sniffles1 and CNVpytor)
filtering SVs | extraction of methylation and phasing info (SVDB | cpgtools)


minimap2, samtools, sniffles 1, BCFtools, CNVpytor, nanopolish, VEP, vcftools, custom database annotation and filtering of variants


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