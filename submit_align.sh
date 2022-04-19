#!/bin/bash -l

#SBATCH -A sens2017106
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J MINIMAP2

module load bioinfo-tools
module load samtools
module load minimap2
module load samtools
module load canu

minimap2 -R '@RG\tID:foo\tSM:bar' -a -t 16 --MD -x map-pb /proj/sens2017106/reference_material/fasta/human_g1k_v37.fasta $1/* | samtools view -Sbh - | samtools sort -m 4G -@16 - > $2
samtools index $2

sbatch sniffles.sh $2
sbatch pytor.sh $2
