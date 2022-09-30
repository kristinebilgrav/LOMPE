# LOMPE
LOng-read Multi-omics PipelinE

Utilize Minimap2, Sniffles, CNVpytor, VEP

Custom analysis of methylation calls from ONT

# RUN
nextflow run main.nf -config < > --folder < > --sample_id < > --style < ont / pb/ pacbio-ccs>

requires;
samtools, canu, minimap2, sniffles (v1), picard

#
