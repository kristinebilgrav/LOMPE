import sys
import statistics

"""
filter methylated calls based on:
log_likelihood score at certain position across reads & above threshold
"""


ratio_dict = {}
methyl_dict = {}
unmethyl_dict = {}
pos_dict = {}
read_dict = {}

chr = ''
filevars= {}
for input in open(sys.argv[1]):
    line= input.rstrip('\n').split('\t')

    if input.startswith('chrom'):
        for field in line:
            filevars[field] = line.index(field)

        continue

    chr = line[filevars['chromosome']] 
    pos_start = line[filevars['start']]
    pos_end = line[filevars['end']]
    read_name = line[filevars['read_name']]
    log_ratio = float(line[filevars['log_lik_ratio']])
    log_methy = float(line[filevars['log_lik_methylated']])
    log_unmeth = float(line[filevars['log_lik_unmethylated']])

    if log_methy > -10 : ##?? more than -10
        continue

    if pos_start not in ratio_dict:
        ratio_dict[pos_start] = []
        methyl_dict[pos_start] = []
        unmethyl_dict[pos_start] = []
        read_dict[pos_start] = []

    if pos_start not in pos_dict:
        pos_dict[pos_start] = pos_end
    else:
        if pos_dict[pos_start] != pos_end:
            print('error', pos_start, pos_end)

    ratio_dict[pos_start].append(log_ratio)
    methyl_dict[pos_start].append(log_methy)
    unmethyl_dict[pos_start].append(log_unmeth)
    read_dict[pos_start].append(read_name)

output = open(sys.argv[2], 'w')
header= ['chr', 'start', 'end', 'log_ratio', 'number_reads', 'reads']
output.write('\t'.join(header) + '\n')
for pos in ratio_dict:
    if len(ratio_dict[pos]) < 2:
        continue
    #if mean log ratio at position above 2, keep
    mean_ratio = sum(ratio_dict[pos])/len(ratio_dict[pos])
    mean_meth =  sum(methyl_dict[pos])/len(methyl_dict[pos])
    mean_unmeth =  sum(unmethyl_dict[pos])/len(unmethyl_dict[pos])
    if mean_ratio > 2:
        #print(mean_ratio,len(ratio_dict[pos]))
        #print(mean_meth)
        #print(mean_unmeth)
        #print(pos)
        #print(pos_dict[pos])
        #print(sorted(ratio_dict[pos]))
        #print(methyl_dict[pos])
        #print(unmethyl_dict[pos])
        myprint = [str(chr), str(pos), str(pos_dict[pos]), str(mean_ratio),  str(len(read_dict[pos])), ';'.join(read_dict[pos])]
        output.write('\t'.join(myprint) + '\n')


