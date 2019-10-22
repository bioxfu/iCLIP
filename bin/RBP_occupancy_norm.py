import sys

gene_sum = {}
gene_site = {}

with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[3]
		cnt = int(lst[4])

		if cnt >= int(sys.argv[2]):
			if ',' not in gene:
				if gene in gene_sum:
					gene_sum[gene] += cnt
					gene_site[gene] += 1
				else:
					gene_sum[gene] = cnt
					gene_site[gene] = 1

with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[3]
		cnt = int(lst[4])

		if cnt >= int(sys.argv[2]):
			if ',' not in gene and gene_site[gene] > 1:
				norm_cnt = cnt / gene_sum[gene] * 100
				print("%s\t%s\t%s\t%s\t%.6f\t%s" % (lst[0], lst[1], lst[2], gene, norm_cnt, lst[5]))
