# LOMPE
LOng-read Multi-omics PipelinE

** Under development **

A pipeline for aligning long read pacbio or nanopore fastq files,
calling SNVs, SVs and methylation, 
with phasing and filtering. 


Workflow: #image
fastq -> bam (minimap2)
snv calling and annotation (bcftools)
phasing and haplotagging of bam file (whatshap)
methylation calling (ont) (nanopolish)
sv calling and cnv calling (sniffles1 and CNVpytor)
filtering SVs | extraction of methylation info (SVDB | nanopolsih/cpgtools)


minimap2, samtools, sniffles 1, BCFtools, CNVpytor, nanopolish, VEP, pb-cpg-tools, custom database annotation and filtering of variants


QC: 
FastQC, picard

# Install

Dependencies: 

Nextflow version 22.10.2

python3 

Docker or singularity

samtools 

bcftools  

FastQ

For ONT methylation calling: hdf5

git clone < repo >

edit config to your needs

# RUN
    nextflow run main.nf -c < config > 
    --input < folder containing 'fastq_pass' folder with fastq files/ samplesheet with folder paths or bam paths (pb) / bam (annotation and methylation calling from pb bam file) > 
    --output < pathto/output/folder> will generate pathto/output/folder/SampleID_out
    --style < ont OR pb>   -with-trace 
    --file < 'fastq' or 'bam' only needed for pb samples>


for ONT samples, a folder called 'fast5_pass' is expected in the same directory as the 'fastq_pass' folder

the samplesheet.csv need to contain a header: SampleID,SamplePath where rows follow with sampleid,folder/orbam/path 
if samplesheet contains PB bam files, a column SampleIndex with bai file path is required

Input:
nanopore: folder containing fastq_pass and fast5_pass folders with fastq.gz / fast5 
pacbio: folder containing fastq_pass with fastq.gz OR bam file

Methylation annotation of bam:
Only performed on ONT samples. 

For SVDB database annotation and filtering:
Local SVDB database with variants to use for filtering

The config file needs to be edited with path to reference genome, containers / paths to program, 
clusteroptions, annotation databases. 
example config file can be found as LOMPE.config

# Output
    ONT: 
    sample_id.phased.haplotagged.bam(.bai) HP and Mm tags 
    sample_id.snv.filtered.vcf
    sample_id.output.VEP.query.filtered.vcf
    sample_id.methylsites.tsv with ref positions, read name and log_likihood of methylation (nanopolish)
    sample_id.methylfrequency.tsv with log_lik_ratio summarized for each position (nanopolish)
    sample_id.wgsmetrics.txt
    sample_id.fastQC.out/

    PB: 
    sample_id.phased.bam(.bai)
    sample_id.snv.filtered.vcf
    sample_id.output.VEP.query.filtered.vcf
    sample_id.wgsmetrics.txt
    sample_id.fastQC.out/

    if bam with methyl sites:
    sample_id.combined/hap1/hap2.denovo.bed pb-cpg-tools for more info
