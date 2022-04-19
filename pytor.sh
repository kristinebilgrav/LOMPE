#!/bin/bash -l
 
#SBATCH -A sens2017106
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J BCF_PYTOR

module load bioinfo-tools bcftools
#bwa index $1

bcftools mpileup -X pacbio-ccs --threads 4 -f /proj/sens2017106/reference_material/fasta/human_g1k_v37.fasta $1 | bcftools call --threads 4 -mv -P0.01 -Ob -o $1.snv.vcf

singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -rd $1
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -gc /proj/sens2017106/reference_material/fasta/human_g1k_v37.fasta
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -his 20000
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -partition 20000
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -call 20000 > $1.pytor.out
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -snp $1.snv.vcf -nofilter
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -baf 20000
singularity exec /proj/sens2017106/nobackup/nano_pore_analysis/cnvpytor_latest.sif cnvpytor -root $1.pytor -call combined 20000 > $1.pytor.combined.combined.out
