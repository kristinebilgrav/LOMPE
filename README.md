# LOMPE
LOng-read Multi-omics PipelinE

* Under development *

A pipeline for aligning long read pacbio or nanopore fastq files,
calling SNVs, SVs and methylation, 
with phasing and filtering.

Returns bam file with haplotags and methylation tags (methylation optional) for visualization and analysis,
vcf files with filtered and annotated snps and svs, 
bed file with positional methylation status (optional)

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

# Install

Dependencies:
Nextflow version 21.10.6
python3
Docker or singularity
samtools 
bcftools  
FastQ
For ONT methylation calling: hdf5

git clone < repo >

edit config to your needs

# RUN
nextflow run main.nf -config < > --fastq_folder < > --fast5_folder < optional > --sample_id < > 
--style < ont OR pb >  
--output <  > -with-trace 

Input:

pacbio (pb) or nanopore (ont) fastq files (gzipped)

Methylation calling:
Only performed on ONT samples. 
Basecalled fast5 files (ont)

For SVDB database annotation and filtering:

Local SVDB database with variants to use for filtering

# Output

returns aligned bam with tags for phasing (HP) and methylation (Mm ; ONT only),
annotated and phased snv file, 
annotated and (given database) filtered SV vcf, 
bed files with phased methylation probabilities (ONT) #implement for PB if present in bam (no align option)

# Installation
git clone

modify config to your needs (executor, reference files, databases)

run!