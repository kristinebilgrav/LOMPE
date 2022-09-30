# LOMPE
LOng-read Multi-omics PipelinE

Utilize Minimap2, Sniffles, CNVpytor, VEP, nanopolish

Custom analysis of methylation calls from ONT

# RUN
nextflow run main.nf -config < > --folder < > --sample_id < > --style < ont OR pb >

requires;
samtools, canu, minimap2, sniffles (v1), picard, pandas?

#
