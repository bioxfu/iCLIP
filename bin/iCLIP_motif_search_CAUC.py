#! /usr/bin/env python
'''
Usage:
    ./iCLIP_motif_search_CAUC.py [aln_reads]
'''

import sys, os

if len(sys.argv) != 2:
    print __doc__
    sys.exit()

reads_file = sys.argv[1]
output_file = 'motif/motifs_CAUC_prop.tsv'
output = open(output_file, 'w')

motifs = {'UCAUC|UUAUC': 0, 'CAUC|UAUC': 0}

total_count = 0
for line in open(reads_file):
    lst = line.strip().split('\t')
    seq, num = lst[1], int(lst[2])
    total_count += num
    if seq.find('TCATC') != -1 or seq.find('TTATC') != -1:
            motifs['UCAUC|UUAUC'] += num 
    if seq.find('CATC') != -1 or seq.find('TATC') != -1:
            motifs['CAUC|UAUC'] += num 

print motifs
print total_count

for motif in motifs:
    output.write('%s\t%.2f%%\n' % (motif, float(motifs[motif])/total_count*100))

output.close()

