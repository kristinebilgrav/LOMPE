# LOMPE
LOng-read Multi-omics PipelinE
nextflow version 21.10.6

Utilize Minimap2, Sniffles, CNVpytor, VEP, nanopolish

Custom analysis of methylation calls from ONT

# RUN
nextflow run main.nf -config < > --fastq_folder < > --fast5_folder < optional > --sample_id < > --style < ont OR pb > -with-trace --output <  >



requires;
pb or ont fastq files
Albacore basecalled fast5 files
samtools, canu, minimap2, sniffles (v1), picard, pandas?, VEP, vcftools, bcftools

#
