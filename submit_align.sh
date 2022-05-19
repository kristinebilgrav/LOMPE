#!/bin/bash -l
#SBATCH -A sens2020021
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J MINIMAP2

module load bioinfo-tools
module load samtools
module load minimap2
module load samtools
module load canu

#pacbio : -x map-pb
#ont: map-ont

minimap2 -R '@RG\tID:foo\tSM:bar' -a -t 16 --MD -x map-ont /sw/data/reference/Homo_sapiens/g1k_v37/downloads/human_g1k_v37.fasta.gz $1/* | samtools view -Sbh - | samtools sort -m 4G -@16 - > $2
samtools index $2

PTH=/proj/sens2020021/ONT/LRpipe
sbatch $PTH/sniffles.sh $2
sbatch $PTH/pytor.sh $2
sbatch $PTH/picard.sh $2
sbatch $PTH/fastqc.sh $2 $2.QC.out
