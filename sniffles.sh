#!/bin/bash -l

#SBATCH -A sens2017106
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J SNIFFLES

#module load bioinfo-tools
#module load samtools
#module load minimap2
#module load samtools
#module load canu


/proj/sens2017106/nobackup/wharf/jesperei/jesperei-sens2017106/Sniffles-master/bin/sniffles-core-1.0.11/sniffles -m $1 -v $1.vcf -l 100 -t 16 -s 3 --genotype --cluster

sbatch annotateVEP37.sh $1.vcf
