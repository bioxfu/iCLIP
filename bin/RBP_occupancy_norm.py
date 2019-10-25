import sys

lib_size = 0
gene_sum = {}
gene_site = {}
gene_len = {}

# gene length
with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[0]
		length = int(lst[1])
		gene_len[gene] = length / 1000

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[3]
		cnt = int(lst[4])

		if cnt >= int(sys.argv[3]):
			if ',' not in gene:
				lib_size += cnt
				if gene in gene_sum:
					gene_sum[gene] += cnt
					gene_site[gene] += 1
				else:
					gene_sum[gene] = cnt
					gene_site[gene] = 1

lib_size = lib_size / 10000000000

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[3]
		cnt = int(lst[4])

		if cnt >= int(sys.argv[3]):
			if ',' not in gene and gene_site[gene] > 1:
				norm_cnt = cnt / gene_sum[gene] / lib_size / gene_len[gene]
				#print("%s\t%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s" % (lst[0], lst[1], lst[2], gene, norm_cnt, lst[5], gene_sum[gene]/gene_site[gene], lib_size, gene_len[gene]))
				print("%s\t%s\t%s\t%s\t%.6f\t%s" % (lst[0], lst[1], lst[2], gene, norm_cnt, lst[5]))
