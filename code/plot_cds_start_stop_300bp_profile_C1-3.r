#! /usr/bin/Rscript

cal_mean_density <- function(plus,minus) {
  plus_mean <- tapply(plus$density, plus$position, mean)
  minus_mean <- rev(tapply(minus$density, minus$position, mean))
  all_mean <- (plus_mean+minus_mean)/2
  return(all_mean)
}

#args <- commandArgs(T)
#input <- args[1]
#output <- args[2]
input <- 'cds_start_stop_300bp_RPKM.10.bed.tag.norm' 
output <- 'cds_start_stop_300bp_RPKM.10.pdf' 

x <- read.table(input,stringsAsFactors=F)
colnames(x) <- c('gene','refseq','strand','end','WT.RPKM','KD.RPKM','position','density')

start_plus  <- x[x$end=='start' & x$strand=='+', ]
start_minus <- x[x$end=='start' & x$strand=='-', ]
stop_plus   <- x[x$end=='stop'  & x$strand=='+', ]
stop_minus  <- x[x$end=='stop'  & x$strand=='-', ]

# start_mean <- cal_mean_density(start_plus, start_minus)
# stop_mean <- cal_mean_density(stop_plus, stop_minus)
# plot(start_mean,type='h',ylab='Normalized iCLIP tag count',xlab="Distance from 5'end of CDS (nt)",xaxt='n')
# axis(1,at=c(0,100,200,300,400,500,600),label=c(-300,-200,-100, 0, 100,200,300))
# plot(stop_mean,type='h',ylab='Normalized iCLIP tag count',xlab="Distance from 3'end of CDS (nt)",xaxt='n')
# axis(1,at=c(0,100,200,300,400,500,600),label=c(-300,-200,-100, 0, 100,200,300))
# barplot(start_mean[291:311],names=c(-10:-1,0,1:10),ylab='Normalized iCLIP tag count',xlab="Distance from 5'end of CDS (nt)")
# barplot(stop_mean[291:311],names=c(-10:-1,0,1:10),ylab='Normalized iCLIP tag count',xlab="Distance from 3'end of CDS (nt)")

start_peak_density <- x$density[x$end=='start' & ((x$strand=='+' & x$position==300) | (x$strand=='-' & x$position==302))] 
start_peak_gene <- x$gene[x$end=='start' & ((x$strand=='+' & x$position==300) | (x$strand=='-' & x$position==302))] 
start_peak_density_rpkm <- x[x$end=='start' & ((x$strand=='+' & x$position==300) | (x$strand=='-' & x$position==302)), c('gene','refseq','density','WT.RPKM','KD.RPKM')] 

# sum(start_peak_density == 0)
# sum(start_peak_density > 0 & start_peak_density <= 0.2)
# sum(start_peak_density > 0.2)

start_peak_gene_C1 <- start_peak_gene[start_peak_density == 0]
start_peak_gene_C2 <- start_peak_gene[start_peak_density > 0 & start_peak_density <= 0.2]
start_peak_gene_C3 <- start_peak_gene[start_peak_density > 0.2]

start_peak_gene_C1.rpkm <- start_peak_density_rpkm[start_peak_density_rpkm$gene %in% start_peak_gene_C1, ]
start_peak_gene_C3.rpkm <- start_peak_density_rpkm[start_peak_density_rpkm$gene %in% start_peak_gene_C3, ]
start_peak_gene_C3.rpkm <- start_peak_gene_C3.rpkm[order(start_peak_gene_C3.rpkm$density,decreasing=T),]
write.table(start_peak_gene_C1.rpkm,'cds_start_peak_gene_C1.rpkm.tsv',row.names=F,quote=F,sep='\t')
write.table(start_peak_gene_C3.rpkm,'cds_start_peak_gene_C3.rpkm.tsv',row.names=F,quote=F,sep='\t')

start_plus_C1 <- start_plus[start_plus$gene %in% start_peak_gene_C1,]
start_plus_C2 <- start_plus[start_plus$gene %in% start_peak_gene_C2,]
start_plus_C3 <- start_plus[start_plus$gene %in% start_peak_gene_C3,]
start_minus_C1 <- start_minus[start_minus$gene %in% start_peak_gene_C1,]
start_minus_C2 <- start_minus[start_minus$gene %in% start_peak_gene_C2,]
start_minus_C3 <- start_minus[start_minus$gene %in% start_peak_gene_C3,]
start_mean_C1 <- cal_mean_density(start_plus_C1, start_minus_C1)
start_mean_C2 <- cal_mean_density(start_plus_C2, start_minus_C2)
start_mean_C3 <- cal_mean_density(start_plus_C3, start_minus_C3)

stop_plus_C1 <- stop_plus[stop_plus$gene %in% start_peak_gene_C1,]
stop_plus_C2 <- stop_plus[stop_plus$gene %in% start_peak_gene_C2,]
stop_plus_C3 <- stop_plus[stop_plus$gene %in% start_peak_gene_C3,]
stop_minus_C1 <- stop_minus[stop_minus$gene %in% start_peak_gene_C1,]
stop_minus_C2 <- stop_minus[stop_minus$gene %in% start_peak_gene_C2,]
stop_minus_C3 <- stop_minus[stop_minus$gene %in% start_peak_gene_C3,]
stop_mean_C1 <- cal_mean_density(stop_plus_C1, stop_minus_C1)
stop_mean_C2 <- cal_mean_density(stop_plus_C2, stop_minus_C2)
stop_mean_C3 <- cal_mean_density(stop_plus_C3, stop_minus_C3)

#pdf(paste0(output,'.pdf'),width=10,heigh=5)
pdf(output,width=10,heigh=5)
par(mfrow=c(1,2))
library(RColorBrewer)
cols <- brewer.pal(3,'Set1')
plot(start_mean_C3,type='h',ylab='iCLIP tag density',xlab="Distance from 5'end of CDS (nt)",xaxt='n',col=cols[1])
axis(1,at=c(0,100,200,300,400,500,600),label=c(-300,-200,-100, 0, 100,200,300))
points(start_mean_C2,type='h',col=cols[2])
points(start_mean_C1,type='h',col=cols[3])
# barplot(start_mean_C3[291:311],names=c(-10:-1,0,1:10),ylab='Normalized iCLIP tag count',xlab="Distance from 5'end of CDS (nt)",col=cols[1],border=NA)
# barplot(start_mean_C2[291:311],yaxt='n',names=NA,add=T,col=cols[2],border=NA)
# barplot(start_mean_C1[291:311],yaxt='n',names=NA,add=T,col=cols[3],border=NA)
bp <- barplot(start_mean_C3[291:310],names=NA,ylab='iCLIP tag density',xlab="",col=cols[1],border=NA)
barplot(start_mean_C2[291:310],yaxt='n',names=NA,add=T,col=cols[2],border=NA)
barplot(start_mean_C1[291:310],yaxt='n',names=NA,add=T,col=cols[3],border=NA)
legend('topright',paste0('C',1:3),fill=rev(cols),bty='n',border=NA)
mtext(c(-10:-1,1:10), at=bp, side=1, cex=0.8)

# 
# #plot(stop_mean_C3,type='h',ylab='Normalized iCLIP tag count',xlab="Distance from 3'end of CDS (nt)",xaxt='n',col=cols[1])
# #axis(1,at=c(0,100,200,300,400,500,600),label=c(-300,-200,-100, 0, 100,200,300))
# plot(stop_mean_C3,type='h',ylab='Normalized iCLIP tag count',xlab="",xaxt='n',col=cols[1])
# axis(1,at=c(0,100,200,300,400,500,600),label=c(-300,-200,-100, 1, 100,200,300))
# points(stop_mean_C2,type='h',col=cols[2])
# points(stop_mean_C1,type='h',col=cols[3])
# # barplot(stop_mean_C3[291:311],names=c(-10:-1,0,1:10),ylab='Normalized iCLIP tag count',xlab="Distance from 3'end of CDS (nt)",col=cols[1],border=NA)
# # barplot(stop_mean_C2[291:311],yaxt='n',names=NA,add=T,col=cols[2],border=NA)
# # barplot(stop_mean_C1[291:311],yaxt='n',names=NA,add=T,col=cols[3],border=NA)
# bp <- barplot(stop_mean_C3[292:311],names=NA,ylab='Normalized iCLIP tag count',xlab="",col=cols[1],border=NA)
# barplot(stop_mean_C2[292:311],yaxt='n',names=NA,add=T,col=cols[2],border=NA)
# barplot(stop_mean_C1[292:311],yaxt='n',names=NA,add=T,col=cols[3],border=NA)
# legend('topright',paste0('C',1:3),fill=rev(cols),bty='n',border=NA)
# mtext(c(-10:-1,1:10), at=bp, side=1)
# 
dev.off()

