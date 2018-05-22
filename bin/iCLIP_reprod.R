#! /usr/bin/env Rscript

args <- commandArgs(T)
input <- args[1]
output <- args[2]
x <- read.table(input, head=T, stringsAsFactors=F)

pdf(output, width=10, heigh=5)
par(mfrow=c(1,2))

for ( f in unique(x$file) ) {
  w0  <- x[x$file==f & x$window==0, 4:6]
  w5  <- x[x$file==f & x$window==5, 4:6] / w0 * 100
  w30 <- x[x$file==f & x$window==30, 4:6] / w0 * 100
  barplot(as.matrix(w5),beside=T,ylim=c(0,100),names=c('c > 0','c > 1','c > 2'),ylab='Percentage of crosslink sites (%)',sub='+/-5bp of crosslink site',main=sub('sites/','',f),cex.axis=1.5,cex.lab=1.5,cex.name=1.5,cex.main=1.5,cex.sub=1.5)
  barplot(as.matrix(w30),beside=T,ylim=c(0,100),names=c('c > 0','c > 1','c > 2'),ylab='Percentage of crosslink sites (%)',sub='+/-30bp of crosslink site',main=sub('sites/','',f),cex.axis=1.5,cex.lab=1.5,cex.name=1.5,cex.main=1.5,cex.sub=1.5)
}

dev.off()

