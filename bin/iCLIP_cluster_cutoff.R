#! /usr/bin/env Rscript

args <- commandArgs(T)
input <- args[1]
output <- args[2]

x <- read.table(input)
tag <- sum(x$V3)
cl <- c(sum(x[,2] > 0), sum(x[,2] > 1), sum(x[,2] > 2), sum(x[,2] > 3), sum(x[,2] > 4))
names(cl) <- c('> 0', '> 1', '> 2', '> 3', '> 4')

pdf(output, width=5, heigh=5)
bp <- barplot(cl/tag,ylab='# clusters / # crosslink sites',xlab='# tags in cluster')
text(bp,cl/tag/2,cl)
dev.off()

