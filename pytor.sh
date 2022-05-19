#!/bin/bash -l
#SBATCH -A sens2020021
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J BCF_PYTOR

module load bioinfo-tools bcftools
#bwa index $1
#pacbio-ccs
#ont

#--annotate FORMAT/AD,FORMAT/DP
#-P0.005

bcftools mpileup -X ont --threads 4 -f /proj/sens2020021/ONT/human_g1k_v37.fasta $1 --annotate FORMAT/AD,FORMAT/DP | bcftools call --threads 4 -mv -P0.01 -Ob -o $1.snv.vcf

singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -rd $1
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -gc /proj/sens2020021/ONT/human_g1k_v37.fasta
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -his 20000 200000
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -partition 20000 200000
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -call 20000 > $1.pytor.out
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -snp $1.snv.vcf -nofilter
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -baf 20000 200000
singularity exec /proj/sens2020021/ONT/cnvpytor_latest.sif cnvpytor -root $1.pytor -call combined 20000 > $1.pytor.combined.out
