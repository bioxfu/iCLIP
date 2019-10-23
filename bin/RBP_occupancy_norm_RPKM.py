import sys

gene_sum = {}
gene_site = {}
gene_rpkm = {}

with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		if lst[0] != 'Name-Gene':
			gene = lst[0]
			rpkm = (float(lst[1]) + float(lst[2])) / 2
			if rpkm > 0:
				gene_rpkm[gene] = rpkm

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[3]
		cnt = int(lst[4])

		if cnt >= int(sys.argv[3]):
			if ',' not in gene:
				if gene in gene_sum:
					gene_sum[gene] += cnt
					gene_site[gene] += 1
				else:
					gene_sum[gene] = cnt
					gene_site[gene] = 1

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[3]
		cnt = int(lst[4])

		if cnt >= int(sys.argv[3]):
			if ',' not in gene and gene_site[gene] > 1 and gene in gene_rpkm:
				norm_cnt = cnt / gene_sum[gene] / gene_rpkm[gene] * 1000
				print("%s\t%s\t%s\t%s\t%.6f\t%s" % (lst[0], lst[1], lst[2], gene, norm_cnt, lst[5]))
				#print("%s\t%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s" % (lst[0], lst[1], lst[2], gene, norm_cnt, lst[5], cnt, gene_sum[gene], gene_rpkm[gene]))
