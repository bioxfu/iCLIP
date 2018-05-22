#! /usr/bin/env python
'''
Usage:
    ./iCLIP_motif_search.py [aln_reads] [motif_list]
'''

import sys, os

if len(sys.argv) != 3:
    print __doc__
    sys.exit()

reads_file, motif_file = sys.argv[1:]
output_file = motif_file +'.prop.tsv'
output = open(output_file, 'w')

motifs = []

for line in open(motif_file):
    motifs.append(line.strip())


motifs_count = [0] * len(motifs)

total_count = 0
for line in open(reads_file):
    lst = line.strip().split('\t')
    seq, num = lst[1], int(lst[2])
    total_count += num
    for i,motif in enumerate(motifs):
        if seq.find(motif) != -1:
            motifs_count[i] += num 

print motifs
print motifs_count
print total_count

for i,motif in enumerate(motifs):
    output.write('%s\t%.2f%%\n' % (motif, float(motifs_count[i])/total_count*100))

output.close()

