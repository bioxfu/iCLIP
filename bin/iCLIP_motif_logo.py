#! /usr/bin/env python
'''
Usage:
    ./iCLIP_motif_logo.py [zscore.tsv] [topNum]
'''

import sys, os

if len(sys.argv) != 3:
    print __doc__
    sys.exit()

input_file, topNum = sys.argv[1:]
fasta_file = os.path.splitext(input_file)[0]+'.fa'
dnd_file = os.path.splitext(input_file)[0]+'.dnd'

seq_lst = []
output = open(fasta_file, 'w')

for line in open(input_file):
    lst = line.strip().split('\t')
    seq_lst.append(lst[0])

for i,seq in enumerate(seq_lst[:20]):
    output.write('>%s\n%s\n' % (i+1, seq))

output.close()

cmd = 'clustalw -INFILE=%s -ALIGN && rm %s %s' % (fasta_file, fasta_file, dnd_file)
os.system(cmd)






