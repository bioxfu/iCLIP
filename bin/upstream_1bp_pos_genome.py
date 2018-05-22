#! /usr/bin/env python

''' the input is bowtie output file
the output is bed file of the m6A site
'''

import sys

for line in sys.stdin:
    line = line.split('\t')
    chrom, strand, start, seq = line[2], line[1], int(line[3]), line[4]
    end = start + len(seq)
    A_start = start - 1
    A_end = A_start + 1
    if strand == '-':
        A_end = end + 1
        A_start = A_end - 1
    id = chrom+':'+str(A_end)+':'+strand
    print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, A_start, A_end, id, seq, strand, len(seq)))

