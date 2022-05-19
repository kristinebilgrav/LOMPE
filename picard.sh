#!/bin/bash -l
 
#SBATCH -A sens2020021
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 4-00:00:00
#SBATCH -J picard

module load bioinfo-tools picard FastQC

java -jar $PICARD_ROOT/picard.jar CollectWgsMetrics --INPUT $1 --OUTPUT $1.wgsmetrics.txt --REFERENCE_SEQUENCE /proj/sens2020021/ONT/human_g1k_v37.fasta --COUNT_UNPAIRED true --MINIMUM_BASE_QUALITY 1 --VALIDATION_STRINGENCY SILENT

