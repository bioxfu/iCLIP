#! /usr/bin/env python

import sys

tot = 0
uniq_aln = 0
mult_aln = 0

with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		if 'Number of input reads' in line:
			tot = int(lst[1])
		elif 'Uniquely mapped reads number' in line:
			uniq_aln = int(lst[1])
		elif 'Number of reads mapped to multiple loci' in line:
			mult_aln += int(lst[1])
		elif 'Number of reads mapped to too many loci' in line:
			mult_aln += int(lst[1])

aln = uniq_aln + mult_aln

print("Number_of_input_reads:\t%s" % tot)
print("Number_of_unique_aligned_reads:\t%s(%.2f%%)" % (uniq_aln, float(uniq_aln) / tot * 100))
print("Number_of_aligned_reads:\t%s(%.2f%%)" % (aln, float(aln) / tot * 100))
