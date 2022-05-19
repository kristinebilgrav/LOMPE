#!/bin/bash -l
#SBATCH -A sens2020021
#SBATCH -p node
#SBATCH -n 6
#SBATCH -t 3-00:00:00
#SBATCH -J SNIFFLES

#module load bioinfo-tools
#module load samtools
#module load minimap2
#module load samtools

module load bioinfo-tools Sniffles


sniffles -m $1 -v $1.vcf -l 100 -t 16 -s 3 --genotype --cluster

PTH=/proj/sens2020021/ONT/LRpipe
sbatch $PTH/annotateVEP37.sh $1.vcf
