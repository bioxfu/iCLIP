#! /usr/bin/env Rscript

args <- commandArgs(T)
feature <- read.table(args[1])
pdf_out <- args[2]

library(RColorBrewer)
pdf(pdf_out)
pie(feature[,2],labels=NA,col=brewer.pal(7,'Set3'),border=F,cex=2,clockwise=T)
legend('topleft',legend=feature[,1],fill=brewer.pal(7,'Set3'),bty='n')
dev.off()


