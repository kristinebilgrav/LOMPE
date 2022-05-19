#!/bin/bash -l
#SBATCH -A sens2020021
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 7-00:00:00
#SBATCH -J VEP

module load bioinfo-tools vep vcftools tabix 

#if not done on files:
#bedtools sort -i file

#cd $TMPDIR

#vcf-sort -c  $1  > $TMPDIR/$2_sorted.vcf
#bgzip $TMPDIR/$2_sorted.vcf
#tabix $TMPDIR/$2_sorted.vcf.gz


#vep --cache --dir $VEP_CACHE --offline  -i $TMPDIR/$2_sorted.vcf.gz -o $TMPDIR/$2.annotated.vcf  --vcf --assembly GRCh37 --per_gene --format vcf --no_stats
vep --cache --dir $VEP_CACHE --offline  -i $1 -o $1.annotated.vcf  --vcf --assembly GRCh37 --per_gene --format vcf --no_stats

#PATH=/proj/sens2017106/kristine/PacBio/
#mv $TMPDIR/$2.annotated.vcf $PATH


