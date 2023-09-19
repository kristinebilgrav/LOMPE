from modbampy import ModBam
import sys

"""
extract methylation from bam files
1: bam file
2: fai
"""

def extract_reads(bam,chrom,start,end):
	for read in bam.reads(chrom,start,end ):
		for pos_mod in read.mod_sites:
			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom,pos_mod.rpos,pos_mod.rpos+1,pos_mod.query_name,read.reference_start,read.reference_end,pos_mod.qpos,pos_mod.strand,pos_mod.qual))
			#break



#ref={"1":1000000}
ref={}
for line in open(sys.argv[2]):
	content=line.strip().split()
	ref[content[0]]=int(content[1])


print("#chromosome\tposition\tposition\tread\tread_start\tread_end\tread_position\tstrand\tquality")
with ModBam(sys.argv[1],) as bam:
	for chrom in ref:
		extract_reads(bam,chrom,1,ref[chrom])	
		#break
