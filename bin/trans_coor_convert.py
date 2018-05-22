#! /usr/bin/env python

import sys

for line in sys.stdin:
	line = line.strip().split('\t')
	dist_to_5p = int(line[2]) - int(line[8])
	seq, chrom, strand = line[4], line[15], line[12]
	if strand == '+':
		A_end = int(line[16]) + dist_to_5p
		A_start = A_end - 1
	else:
		A_start = int(line[17]) - dist_to_5p
		A_end = A_start + 1
	id = chrom+':'+str(A_end)+':'+strand
	print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrom, A_start, A_end, id, seq, strand, len(seq)))

