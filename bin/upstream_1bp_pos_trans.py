#! /usr/bin/env python

''' the input is bowtie output file
the output is bed file of the m6A site
'''

import sys

for line in sys.stdin:
    line = line.split('\t')
    trans, strand, start, seq = line[2], line[1], int(line[3]), line[4]
    A_start = start - 1
    A_end = A_start + 1
    if strand == '+' and A_start > 0:
        id = trans+':'+str(A_end)+':'+strand
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (trans, A_start, A_end, id, seq, strand, len(seq)))

