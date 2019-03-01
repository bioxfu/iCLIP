import sys
import re

flank = int(sys.argv[2])

with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		chrom, start, end, strand = re.split('[:\-\(\)]', lst[0])[:4]
		start = int(start)
		end = int(end)
		seq = lst[1]
		seq_len = len(seq)
		motif = 'GCATG|TGCAT'

		for m in re.finditer(motif, seq):
			sp = m.span()
			if strand == '+':
				motif_start = start + sp[0] - flank
				motif_end = start + sp[1] + flank
			else:
				motif_start = start + (seq_len - sp[1]) - flank
				motif_end = start + (seq_len - sp[0]) + flank
			if motif_start > 0:
				print(chrom+"\t"+str(motif_start)+"\t"+str(motif_end)+"\t"+lst[0])

