#! /usr/bin/env Rscript
library(RColorBrewer)
col_set <- brewer.pal(6, 'Set3')

args <- commandArgs(T)
input1 <- args[1]
input2 <- args[2]
output <- args[3]
# input1 <- './merge/merge_ALL_crosslink_sites_sig_summary.tab'
# input2 <- './merge/merge_ALL_crosslink_sites_summary.tab'
# output <- './figure/crosslink_sites_distr.pdf'

x <- read.table(input1, row.names = 1, header = T)
y <- read.table(input2, row.names = 1, header = T)
n <- c('CDS','UTR3','UTR5','intron','intergenic','ncRNA')

pdf(output, width=8, heigh=8)
par(mfrow=c(2,2))
pie(x[n, 'sites'], lab = paste0(n, '(', round(x[n, 'sites']/sum(x[n, 'sites'])*100, 1), '%)'), col=col_set, border = 'white', main = 'sites (FDR < 0.05)')
pie(x[n, 'reads'], lab = paste0(n, '(', round(x[n, 'reads']/sum(x[n, 'reads'])*100, 1), '%)'), col=col_set, border = 'white', main='reads (FDR < 0.05)')
pie(y[n, 'sites'], lab = paste0(n, '(', round(y[n, 'sites']/sum(y[n, 'sites'])*100, 1), '%)'), col=col_set, border = 'white', main='sites (all)')
pie(y[n, 'reads'], lab = paste0(n, '(', round(y[n, 'reads']/sum(y[n, 'reads'])*100, 1), '%)'), col=col_set, border = 'white', main='reads (all)')
dev.off()
