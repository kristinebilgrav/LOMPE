#!/bin/bash -l

#SBATCH -A sens2017106
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 4-00:00:00
#SBATCH -J FastQC

module load bioinfo-tools FastQC

zcat $1 | head -5000000 > $2.5m.fastq
mkdir $2.5m.fastQC.out
fastqc --threads 16 -o $2.5m.fastQC.out $2.5m.fastq


#1: one fastq
#2: output name
