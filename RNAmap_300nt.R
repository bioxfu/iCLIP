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
output <- paste0(prefix, "RNAmap_filt_300nt.pdf")

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

bin300 <- c(grep('^1.', names(si), value = T)[1:35],
            rev(rev(grep('^2.', names(si), value = T))[1:35]),
            grep('^3.', names(si), value = T)[1:35],
            rev(rev(grep('^4.', names(si), value = T))[1:35]))
en <- en[bin300]
si <- si[bin300]
bg <- bg[bin300]
si_pvalue <- si_pvalue[bin300]
en_pvalue <- en_pvalue[bin300]

pdf(output,width=10,heigh=5)
y_min <- min(c(en, si, bg))
y_max <- max(c(en, si, bg))
# y_range = c(y_min-10, y_max+20)
y_range = c(y_min, y_max*1.3)
#plot(-1000,-1000,xlab='',ylab='Occurrence (%)',ylim=y_range,xlim=c(0,332),bty='n',xaxt='n',yaxt='n')
plot(-1000,-1000,xlab='',ylab='',ylim=y_range,xlim=c(0,155),bty='n',xaxt='n',yaxt='n')
mtext('Occurrence (%)',side=2,at=0,line=2)
# axis(2,at=pretty(c(0,y_max)),label=pretty(c(0,y_max)))
# axis(2,at=pretty(c(y_min,0)),label=-1*pretty(c(y_min,0)))
axis(2,at=c(-40,-20,0,20,40),label=c(40,20,0,20,40))
polygon(c(1,1:35,35), c(0,en[1:35],0),col=red,border=NA)
polygon(c(41,41:75,75), c(0,en[36:70],0),col=red,border=NA)
polygon(c(81,81:115,115), c(0,en[71:105],0),col=red,border=NA)
polygon(c(121,121:155,155), c(0,en[106:140],0),col=red,border=NA)
polygon(c(1,1:35,35), c(0,si[1:35],0),col=blue,border=NA)
polygon(c(41,41:75,75), c(0,si[36:70],0),col=blue,border=NA)
polygon(c(81,81:115,115), c(0,si[71:105],0),col=blue,border=NA)
polygon(c(121,121:155,155), c(0,si[106:140],0),col=blue,border=NA)

polygon(c(1,1:35,35), c(0,bg[1:35],0),col=gray_col,border=NA)
polygon(c(41,41:75,75), c(0,bg[36:70],0),col=gray_col,border=NA)
polygon(c(81,81:115,115), c(0,bg[71:105],0),col=gray_col,border=NA)
polygon(c(121,121:155,155), c(0,bg[106:140],0),col=gray_col,border=NA)

polygon(c(1,1:35,35), c(0,-bg[1:35],0),col=gray_col,border=NA)
polygon(c(41,41:75,75), c(0,-bg[36:70],0),col=gray_col,border=NA)
polygon(c(81,81:115,115), c(0,-bg[71:105],0),col=gray_col,border=NA)
polygon(c(121,121:155,155), c(0,-bg[106:140],0),col=gray_col,border=NA)

points(1:35, si_pvalue[1:35],pch=16,cex=1)
points(41:75, si_pvalue[36:70],pch=16,cex=1)
points(81:115, si_pvalue[71:105],pch=16,cex=1)
points(121:155, si_pvalue[106:140],pch=16,cex=1)

points(1:35, en_pvalue[1:35],pch=16,cex=1)
points(41:75, en_pvalue[36:70],pch=16,cex=1)
points(81:115, en_pvalue[71:105],pch=16,cex=1)
points(121:155, en_pvalue[106:140],pch=16,cex=1)

segments(1+4, y_min, 1+4, y_max*0.8, col=green)
segments(75-4, y_min, 75-4, y_max, col=green)
segments(81+4, y_min, 81+4, y_max, col=green)
segments(155-4, y_min, 155-4, y_max, col=green)
segments(1+4, y_max*1.3, 155-4, y_max*1.3, lwd=2)

abline(h=0)
points(c(seq(5,35,10), seq(41,71,10), seq(85,115,10), seq(121,151,10)), rep(0,16), pch=3)
segments(35, 0, 41, 0, lty = 2, col = 'white')
segments(75, 0, 81, 0, lty = 2, col = 'white')
segments(115, 0, 121, 0, lty = 2, col = 'white')


segments(5, y_max*0.95, 15, y_max*0.95)
text(10, y_max*1.01, '100 nt')
polygon(c(75-4,75-4,81+4), c(y_max*1.25, y_max*1.35, y_max*1.25), col=red, border=NA)
polygon(c(75-4,81+4,81+4), c(y_max*1.35, y_max*1.35, y_max*1.25), col=blue, border=NA)
rect(1+4-10, y_max*1.25, 1+4, y_max*1.35, col='gray', border=NA)
rect(155-4, y_max*1.25, 155-4+10, y_max*1.35, col='gray', border=NA)

text(80,y_max,paste0('Enhanced alternative exons (n = ',en_exon,')'),pos=3,col=red)
mtext(c("5' SS","3' SS","5' SS","3' SS"),side=1,at=c(1+3,165-3,171+3,335-3),adj=c(0,1,0,1))
mtext(paste0('Silenced alternative exons (n = ',si_exon,')'),side=1,col=blue,line=1)
mtext(paste0('Unaffected cassette exons (n = ',bg_exon,')'),side=1,col=gray(0.2),line=2)

dev.off()

