#! /usr/bin/Rscript

library(RColorBrewer)
cols <- brewer.pal(3, 'Set1')
red <- cols[1]
blue <- cols[2]
green <- cols[3]
#gray_col <- gray(1, alpha=0.5)
gray_col <- rgb(190,190,190,120,max=255)

init <- function(file) {
  d <- read.table(file,row.names=1)
  x <- rep(0,320)
  names(x)  <- c(paste(1, seq(-3,76),sep='.'), paste(2, seq(-75,4),sep='.'), paste(3, seq(-3,76),sep='.'), paste(4, seq(-75,4),sep='.'))
  x[names(x) %in% rownames(d)] = d$V2
  return(x)
}

args <- commandArgs(T)
prefix <- args[1]
en_exon <- as.numeric(args[2])
si_exon <- as.numeric(args[3])
bg_exon <- as.numeric(args[4])
p_cutoff <- as.numeric(args[5])

# prefix <- "/data1/HJY/Project/QKI/RNAmap/data"
# en_exon <- as.numeric(244)
# si_exon <- as.numeric(207)
# bg_exon <- as.numeric(9740)
# p_cutoff <- as.numeric(0.000001)

en_f <- init(paste0(prefix, "splice_sig_rnamap_en_f.filt"))
en_r <- init(paste0(prefix, "splice_sig_rnamap_en_r.filt"))
si_f <- init(paste0(prefix, "splice_sig_rnamap_si_f.filt"))
si_r <- init(paste0(prefix, "splice_sig_rnamap_si_r.filt"))
bg_f <- init(paste0(prefix, "splice_nonsig_rnamap_f.filt"))
bg_r <- init(paste0(prefix, "splice_nonsig_rnamap_r.filt"))
output <- paste0(prefix, "RNAmap_filt_800nt.pdf")

en <- (en_f + rev(en_r))/en_exon * 100
si <- (si_f + rev(si_r))/si_exon * 100 * -1
bg <- (bg_f + rev(bg_r))/bg_exon * 100

en_2 <- cbind(en_f + rev(en_r), en_exon, bg_f + rev(bg_r), bg_exon)
en_pvalue = apply(en_2, 1, function(x) {fisher.test(matrix(x,nrow=2,ncol=2))$p.value})
en_sig = en_pvalue < p_cutoff
en_pvalue[en_sig] = en[en_sig]
en_pvalue[!en_sig] = 100000
en_pvalue[abs(en_pvalue)<10] = 100000
print(sum(en_pvalue < 100000))

si_2 <- cbind(si_f + rev(si_r), si_exon, bg_f + rev(bg_r), bg_exon)
si_pvalue = apply(si_2, 1, function(x) {fisher.test(matrix(x,nrow=2,ncol=2))$p.value})
si_sig = si_pvalue < p_cutoff
si_pvalue[si_sig] = si[si_sig]
si_pvalue[!si_sig] = 100000
si_pvalue[abs(si_pvalue)<10] = 100000
print(sum(si_pvalue < 100000))

pdf(output,width=10,heigh=5)
y_min <- min(c(en, si, bg))
y_max <- max(c(en, si, bg))
# y_range = c(y_min-10, y_max+20)
y_range = c(y_min, y_max+20)
#plot(-1000,-1000,xlab='',ylab='Occurrence (%)',ylim=y_range,xlim=c(0,332),bty='n',xaxt='n',yaxt='n')
plot(-1000,-1000,xlab='',ylab='',ylim=y_range,xlim=c(0,332),bty='n',xaxt='n',yaxt='n')
mtext('Occurrence (%)',side=2,at=0,line=2)
# axis(2,at=pretty(c(0,y_max)),label=pretty(c(0,y_max)))
# axis(2,at=pretty(c(y_min,0)),label=-1*pretty(c(y_min,0)))
axis(2,at=c(-20,-10,0,10,20),label=c(20,10,0,10,20))
polygon(c(1,1:80,80), c(0,en[1:80],0),col=red,border=NA)
polygon(c(86,86:165,165), c(0,en[81:160],0),col=red,border=NA)
polygon(c(171,171:250,250), c(0,en[161:240],0),col=red,border=NA)
polygon(c(256,256:335,335), c(0,en[241:320],0),col=red,border=NA)
polygon(c(1,1:80,80), c(0,si[1:80],0),col=blue,border=NA)
polygon(c(86,86:165,165), c(0,si[81:160],0),col=blue,border=NA)
polygon(c(171,171:250,250), c(0,si[161:240],0),col=blue,border=NA)
polygon(c(256,256:335,335), c(0,si[241:320],0),col=blue,border=NA)

polygon(c(1,1:80,80), c(0,bg[1:80],0),col=gray_col,border=NA)
polygon(c(86,86:165,165), c(0,bg[81:160],0),col=gray_col,border=NA)
polygon(c(171,171:250,250), c(0,bg[161:240],0),col=gray_col,border=NA)
polygon(c(256,256:335,335), c(0,bg[241:320],0),col=gray_col,border=NA)
polygon(c(1,1:80,80), c(0,-bg[1:80],0),col=gray_col,border=NA)
polygon(c(86,86:165,165), c(0,-bg[81:160],0),col=gray_col,border=NA)
polygon(c(171,171:250,250), c(0,-bg[161:240],0),col=gray_col,border=NA)
polygon(c(256,256:335,335), c(0,-bg[241:320],0),col=gray_col,border=NA)

points(1:80, si_pvalue[1:80],pch=16,cex=1)
points(86:165, si_pvalue[81:160],pch=16,cex=1)
points(171:250, si_pvalue[161:240],pch=16,cex=1)
points(256:335, si_pvalue[241:320],pch=16,cex=1)
points(1:80, en_pvalue[1:80],pch=16,cex=1)
points(86:165, en_pvalue[81:160],pch=16,cex=1)
points(171:250, en_pvalue[161:240],pch=16,cex=1)
points(256:335, en_pvalue[241:320],pch=16,cex=1)

segments(1+4,y_min,1+4,y_max*0.6, col=green)
segments(165-4,y_min,165-4,y_max, col=green)
segments(171+4,y_min,171+4,y_max, col=green)
segments(335-4,y_min,335-4,y_max, col=green)
segments(1+4,y_max+8,335-4,y_max+8, lwd=2)

abline(h=0)
points(c(seq(1,80,20), 80, seq(86,165,20), 165, seq(171,250,20), 250, seq(256,335,20), 335),rep(0,20),pch=3)
segments(1,y_max-5,21,y_max-5)
arrows(c(1,5),y_max-6,c(1,5),y_max-4,length=0)
arrows(c(5,9),y_max-6,c(5,9),y_max-4,length=0)
arrows(c(9,13),y_max-6,c(9,13),y_max-4,length=0)
arrows(c(13,17),y_max-6,c(13,17),y_max-4,length=0)
arrows(c(17,21),y_max-6,c(17,21),y_max-4,length=0)
text(12,y_max,'200 nt')
polygon(c(165-4,165-4,171+4),c(y_max+8-2,y_max+8+2,y_max+8-2),col=blue,border=NA)
polygon(c(165-4,171+4,171+4),c(y_max+8+2,y_max+8+2,y_max+8-2),col=red,border=NA)
rect(1+4-10,y_max+8-2,1+4,y_max+8+2,col='gray',border=NA)
rect(335-4,y_max+8-2,335-4+10,y_max+8+2,col='gray',border=NA)

text(165+3,y_max,paste0('Enhanced alternative exons (n = ',en_exon,')'),pos=3,col=red)
mtext(c("5' SS","3' SS","5' SS","3' SS"),side=1,at=c(1+3,165-3,171+3,335-3),adj=c(0,1,0,1))
mtext(paste0('Silenced alternative exons (n = ',si_exon,')'),side=1,col=blue,line=1)
mtext(paste0('Unaffected cassette exons (n = ',bg_exon,')'),side=1,col=gray(0.2),line=2)

dev.off()

