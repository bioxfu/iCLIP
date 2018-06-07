#! /usr/bin/Rscript

library(RColorBrewer)
cols <- brewer.pal(4,'Set1')

input1 <- 'refGene.rnaseq.RPKM.10.5UTR.50bins.bed.tag.norm.summary'
input2 <- 'refGene.rnaseq.RPKM.10.CDS.100bins.bed.tag.norm.summary'
input3 <- 'refGene.rnaseq.RPKM.10.3UTR.50bins.bed.tag.norm.summary'
output <- 'mRNA_iCLIP_profile.pdf'

x1 <- read.table(input1,stringsAsFactors=F)
x2 <- read.table(input2,stringsAsFactors=F)
x3 <- read.table(input3,stringsAsFactors=F)
v <- c(x1$V2, x2$V2, x3$V2)  

down <- min(v)
up <- max(v)
heigh <- up - down

pdf(output,height=5,width=5)
plot(v,type='h',col=cols[1],lwd=3,xlab='',ylab='',xaxt='n',ylim=c(down-2*heigh/10, up),bty='n',yaxt='n')
axis(2,at=seq(0,up,0.01),label=seq(0,up,0.01))
rect(0, down-3.5*heigh/20, 50, down-2.5*heigh/20,col=cols[2],border=NA)
rect(50, down-3.5*heigh/20, 150, down-2.5*heigh/20,col=cols[3],border=NA)
rect(150, down-3.5*heigh/20, 200, down-2.5*heigh/20,col=cols[4],border=NA)
abline(h=0)
mtext('crosslink site density',side=2,line=2.2)
points(seq(0,300,25),rep(down-4*heigh/20,13),pch='.',cex=3)
mtext(c(0,50),side=1,at=c(0,25),las=1,col=cols[2])
mtext(seq(0,75,25),side=1,at=seq(0,75,25)+50,las=1,col=cols[3])
mtext(seq(0,100,50),side=1,at=seq(0,50,25)+150,las=1,col=cols[4])
mtext('% of transcribed region',side=1,line=1.5)
legend('topleft',c("5'UTR","CDS","3'UTR"),fill=cols[2:4],bty='n',border=NA)
dev.off()

