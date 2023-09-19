import sys

"""
filters LOMPE output based on frequency, quality metrics (PASS, ) ,
VEP annotations, 
"""
output = open(sys.argv[2], 'w')

for line in open(sys.argv[1]):

    if line.startswith('#'):
        header = 'pass'
        output.write(line)
        continue    

    line = line.rstrip('\n')

    # SV length
    
    svlen = abs(int(line.split('SVLEN=')[-1].split(';')[0]))
    if svlen > 1000:
        length = 'pass'
    else:
        length = 'fail'

    if 'BND' in line:
        length = 'pass'

    # FRQ
    if 'FRQ' in line:
        af = float(line.split('FRQ=')[-1].split('\t')[0])
        if af < 0.05:
            frq = 'pass'
        else:
            frq ='fail'
            continue     
    else:
        frq = 'pass'


    if length == 'pass' and frq == 'pass':
        output.write(line + '\n')

    else:
        continue
