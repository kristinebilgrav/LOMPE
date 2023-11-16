
![lompe_logo](https://user-images.githubusercontent.com/77359122/231418667-484df799-c9a6-4e0c-a691-872704068608.png)

# LOng-read Multi-omics PipelinE


A pipeline for aligning long read pacbio or nanopore fastq files,
calling SNVs, SVs and methylation, 
with phasing and filtering. 

QC: 
- Picard
- fastQC

Alignment:
- minimap2

SNV calling:
- deepvariant 
- bcftools (for lower quality data)

SV calling:
- Sniffles v1

CNV calling:
- CNV pythor

Expansion repeats:
- TRGT (https://github.com/PacificBiosciences/trgt)

Phasing:
- WhatsHap

Methylation:
- PB: pb-cpg-tools
- ONT: specify pore; for < R10 nanopolish, R10


Workflow: #image / clean
fastq -> bam (minimap2)
snv calling and annotation (bcftools)
phasing and haplotagging of bam file (whatshap)
methylation calling (ont) (nanopolish)
filtering SVs | extraction of methylation info (SVDB | nanopolsih/cpgtools)



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
    nextflow run kristinebilgrav/lompe 
    
    REQUIRED:
    -c < config > 
    --input < ONT: folder containing 'fastq_pass' with fastq files/ samplesheet with folder paths PB: samplesheet with paths to bam or fastq file  > 
    --output < pathto/output/folder> will generate pathto/output/folder/SampleID_out
    --style < ont OR pb>   
    --file < 'fastq' OR 'bam' OR 'ubam'>

    OPTIONAL:
    -with-trace < creates log file >
    --bcftools < if for some reason want to use bcftools instead of deepvariant >



the samplesheet.csv need to contain a header: SampleID,SamplePath where rows follow with sampleid,folder/orbam/path 

for PB fastq files a path to the fastq file (gz) is expected, OR a bam file with methylation tags incorported

for ONT samples, a folder called 'fast5_pass' is expected in the same directory as the 'fastq_pass' folder OR path to folder with ubams to remap to one

Input:
nanopore: folder containing fastq_pass and fast5_pass folders with fastq.gz / fast5  OR path to folder containing ubams 
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
    
    

    
