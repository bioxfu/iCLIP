#! /usr/bin/Rscript

library(RColorBrewer)
cols <- brewer.pal(4,'Set1')

x <- read.table('lincRNA.rnaseq.RPKM.1.100bins.bed.tag.norm',stringsAsFactors=F)
y <- tapply(x$V8, x$V1,mean)
pdf('lincRNA_iCLIP_outlier.pdf',height=5,width=5)
boxplot(y,ylab='Average normalized iCLIP tag count',main=paste0('outlier:',sub('\\..+','',names(y)[y>1])))
dev.off()

v <- tapply(x$V8, x$V6,mean)  
filt <- names(y)[y<1]
x_filt <- x[x$V1 %in% filt,]
v_filt <- tapply(x_filt$V8, x_filt$V6,mean)  

pdf('lincRNA_iCLIP_profile_with_outlier.pdf',height=5,width=5)
par(mar=c(0,4,1,2))
layout(matrix(1:2,nrow=2,ncol=1),height = c(5,1))
plot(v,type='h',col=cols[1],lwd=3,xlab='',ylab='crosslink site density',xaxt='n',ylim=c(0, 0.25),bty='n')
legend('topleft','lincRNA',fill=cols[2],bty='n',border=NA)
par(mar=c(0,4,0,2))
plot(-10,-10,xlim=c(0,1),ylim=c(0,0.1),ylab='',xlab='',xaxt='n',yaxt='n',bty='n')
rect(0, 0.09, 1, 1,col=cols[2],border=NA)
points(seq(0,1,0.25),rep(0.08,5),pch='.',cex=3)
text(seq(0,1,0.25),rep(0.06,5),seq(0,100,25),las=1,col=cols[2])
text(0.5,0.03,'% of transcribed region')
dev.off()


pdf('lincRNA_iCLIP_profile_without_outlier.pdf',height=5,width=5)
par(mar=c(0,4,1,2))
layout(matrix(1:2,nrow=2,ncol=1),height = c(5,1))
plot(v_filt,type='h',col=cols[1],lwd=3,xlab='',ylab='crosslink site density',xaxt='n',ylim=c(0, 0.04),bty='n',yaxt='n')
axis(2,seq(0,0.04,0.01),seq(0,0.04,0.01))
legend('topleft','lincRNA',fill=cols[2],bty='n',border=NA)
par(mar=c(0,4,0,2))
plot(-10,-10,xlim=c(0,1),ylim=c(0,0.1),ylab='',xlab='',xaxt='n',yaxt='n',bty='n')
rect(0, 0.09, 1, 1,col=cols[2],border=NA)
points(seq(0,1,0.25),rep(0.08,5),pch='.',cex=3)
text(seq(0,1,0.25),rep(0.06,5),seq(0,100,25),las=1,col=cols[2])
text(0.5,0.03,'% of transcribed region')
dev.off()
