import sys
import statistics

"""
filter methylated calls based on:
log_likelihood score at certain position across reads & above threshold
"""


ratio_dict = {}
methyl_dict = {}
unmethyl_dict = {}

chr = ''
for input in open(sys.argv[1]):
    line= input.rstrip('\n').split('\t')

    if input.startswith('chrom'):
        chr_inx = line.index('chromosome')
        read_name_inx = line.index('read_name')

        start_inx = line.index('start')
        end_inx = line.index('end')

        log_ratio_inx = line.index('log_lik_ratio')
        log_methy_inx = line.index('log_lik_methylated')
        log_unmeth_inx = line.index('log_lik_unmethylated')
        continue

    chr = line[chr_inx] 
    pos_start = line[start_inx]
    pos_end = line[end_inx]
    read_name = line[read_name_inx]
    log_ratio = float(line[log_ratio_inx])
    log_methy = float(line[log_methy_inx])
    log_unmeth = float(line[log_unmeth_inx])

    if pos_start not in ratio_dict:
        ratio_dict[pos_start] = []
        methyl_dict[pos_start] = []
        unmethyl_dict[pos_start] = []
    
    ratio_dict[pos_start].append(log_ratio)
    methyl_dict[pos_start].append(log_methy)
    unmethyl_dict[pos_start].append(log_unmeth)


for pos in ratio_dict:
    mean = sum(ratio_dict[pos])/len(ratio_dict[pos])
    if mean > 10:
        print(mean)
        print(pos)
        print(sorted(ratio_dict[pos]))
        print(methyl_dict[pos])
        print(unmethyl_dict[pos])
        quit()
