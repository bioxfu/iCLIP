#! /usr/bin/env python

import sys

rid = {}
uniq_aln = 0
un_aln = 0

for line in sys.stdin:
	lst = line.strip().split('\t')
	rid[lst[0]] = 1
	if lst[11] == 'NH:i:0':
		un_aln += 1
	elif lst[11] == 'NH:i:1':
		uniq_aln += 1

tot = len(rid.keys())
aln = tot - un_aln

print("Number_of_input_reads:\t%s" % tot)
print("Number_of_unique_aligned_reads:\t%s(%.2f%%)" % (uniq_aln, float(uniq_aln) / tot * 100))
print("Number_of_aligned_reads:\t%s(%.2f%%)" % (aln, float(aln) / tot * 100))
