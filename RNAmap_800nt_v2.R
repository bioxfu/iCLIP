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

# prefix <- "RNAmap/data/"
# en_exon <- as.numeric(244)
# si_exon <- as.numeric(207)

get_value <- function(x) {
  en_f <- init(paste0(prefix, x, "_splice_sig_rnamap_en_f.filt"))
  en_r <- init(paste0(prefix, x, "_splice_sig_rnamap_en_r.filt"))
  si_f <- init(paste0(prefix, x, "_splice_sig_rnamap_si_f.filt"))
  si_r <- init(paste0(prefix, x, "_splice_sig_rnamap_si_r.filt"))
  en <- en_f + rev(en_r)
  si <- si_f + rev(si_r)
  return(list(en=en, si=si))
}

wt1 <- get_value('WT1')
wt2 <- get_value('WT2')
wt3 <- get_value('WT3')
ko1 <- get_value('KO1')
ko2 <- get_value('KO2')
ko3 <- get_value('KO3')

wt_en <- (wt1$en + wt2$en + wt3$en ) / 3
ko_en <- (ko1$en + ko2$en + ko3$en ) / 3
wt_si <- (wt1$si + wt2$si + wt3$si ) / 3
ko_si <- (ko1$si + ko2$si + ko3$si ) / 3


# dfm_en <- as.data.frame(cbind(wt1$en, wt2$en, wt3$en, ko1$en, ko2$en, ko3$en))
# dfm_si <- as.data.frame(cbind(wt1$si, wt2$si, wt3$si, ko1$si, ko2$si, ko3$si))
# dfm_en$p <- apply(dfm_en, 1, function(x){t.test(x[1:3], x[4:6], alternative = 'greater')$p.value})
# dfm_si$p <- apply(dfm_si, 1, function(x){t.test(x[1:3], x[4:6], alternative = 'greater')$p.value})
# dfm_en$p[is.nan(dfm_en$p)] <- 1
# dfm_si$p[is.nan(dfm_si$p)] <- 1


pdf(paste0(prefix, "RNAmap_800nt_enhanced.pdf"), width=10,heigh=5)
y_min <- min(c(wt_en, ko_en))
y_max <- max(c(wt_en, ko_en))
y_range = c(0, y_max*1.3)
plot(-1000, -1000, xlab='', ylab='Total normalized counts in each bin', ylim=y_range, xlim=c(0,332), bty='n', xaxt='n')
lines(c(1:80), wt_en[1:80], col=red)
lines(c(86:165), wt_en[81:160], col=red)
lines(c(171:250), wt_en[161:240], col=red)
lines(c(256:335), wt_en[241:320], col=red)
lines(c(1:80), ko_en[1:80], col=blue)
lines(c(86:165), ko_en[81:160], col=blue)
lines(c(171:250), ko_en[161:240], col=blue)
lines(c(256:335), ko_en[241:320], col=blue)

segments(1+4, y_min, 1+4, y_max*0.8, col=green)
segments(165-4, y_min, 165-4, y_max, col=green)
segments(171+4, y_min, 171+4, y_max, col=green)
segments(335-4, y_min, 335-4, y_max, col=green)
segments(1+4, y_max*1.1, 335-4, y_max*1.1, lwd=2)

abline(h=0)
points(c(seq(1,80,20), 80, seq(86,165,20), 165, seq(171,250,20), 250, seq(256,335,20), 335),rep(0,20),pch=3)
segments(1, y_max*0.95, 21, y_max*0.95)
arrows(c(1,5), y_max*0.93, c(1,5), y_max*0.97, length=0)
arrows(c(5,9), y_max*0.93, c(5,9), y_max*0.97, length=0)
arrows(c(9,13), y_max*0.93, c(9,13), y_max*0.97, length=0)
arrows(c(13,17), y_max*0.93, c(13,17), y_max*0.97, length=0)
arrows(c(17,21), y_max*0.93, c(17,21), y_max*0.97, length=0)
text(12, y_max*1.01, '200 nt')
polygon(c(165-4,165-4,171+4), c(y_max*1.05, y_max*1.15, y_max*1.05), col=blue, border=NA)
polygon(c(165-4,171+4,171+4), c(y_max*1.15, y_max*1.15, y_max*1.05), col=red, border=NA)
rect(1+4-10, y_max*1.05, 1+4, y_max*1.15, col='gray', border=NA)
rect(335-4, y_max*1.05, 335-4+10, y_max*1.15, col='gray', border=NA)

mtext(paste0('Enhanced alternative exons (n = ',en_exon,')'), side = 3)
mtext(c("5' SS","3' SS","5' SS","3' SS"), side=1, at=c(1+3, 165-3, 171+3, 335-3), adj=c(0,1,0,1))

legend('topright', c('WT', 'KO'), col=c(red, blue), bty='n',  lwd=1.5, ncol = 2)
dev.off()


pdf(paste0(prefix, "RNAmap_800nt_silenced.pdf"), width=10,heigh=5)
y_min <- min(c(wt_si, ko_si))
y_max <- max(c(wt_si, ko_si))
y_range = c(0, y_max*1.3)
plot(-1000, -1000, xlab='', ylab='Total normalized counts in each bin', ylim=y_range, xlim=c(0,332), bty='n', xaxt='n')
lines(c(1:80), wt_si[1:80], col=red)
lines(c(86:165), wt_si[81:160], col=red)
lines(c(171:250), wt_si[161:240], col=red)
lines(c(256:335), wt_si[241:320], col=red)
lines(c(1:80), ko_si[1:80], col=blue)
lines(c(86:165), ko_si[81:160], col=blue)
lines(c(171:250), ko_si[161:240], col=blue)
lines(c(256:335), ko_si[241:320], col=blue)

segments(1+4, y_min, 1+4, y_max*0.8, col=green)
segments(165-4, y_min, 165-4, y_max, col=green)
segments(171+4, y_min, 171+4, y_max, col=green)
segments(335-4, y_min, 335-4, y_max, col=green)
segments(1+4, y_max*1.1, 335-4, y_max*1.1, lwd=2)

abline(h=0)
points(c(seq(1,80,20), 80, seq(86,165,20), 165, seq(171,250,20), 250, seq(256,335,20), 335),rep(0,20),pch=3)
segments(1, y_max*0.95, 21, y_max*0.95)
arrows(c(1,5), y_max*0.93, c(1,5), y_max*0.97, length=0)
arrows(c(5,9), y_max*0.93, c(5,9), y_max*0.97, length=0)
arrows(c(9,13), y_max*0.93, c(9,13), y_max*0.97, length=0)
arrows(c(13,17), y_max*0.93, c(13,17), y_max*0.97, length=0)
arrows(c(17,21), y_max*0.93, c(17,21), y_max*0.97, length=0)
text(12, y_max*1.01, '200 nt')
polygon(c(165-4,165-4,171+4), c(y_max*1.05, y_max*1.15, y_max*1.05), col=blue, border=NA)
polygon(c(165-4,171+4,171+4), c(y_max*1.15, y_max*1.15, y_max*1.05), col=red, border=NA)
rect(1+4-10, y_max*1.05, 1+4, y_max*1.15, col='gray', border=NA)
rect(335-4, y_max*1.05, 335-4+10, y_max*1.15, col='gray', border=NA)

mtext(paste0('Silenced alternative exons (n = ',si_exon,')'), side = 3)
mtext(c("5' SS","3' SS","5' SS","3' SS"), side=1, at=c(1+3, 165-3, 171+3, 335-3), adj=c(0,1,0,1))

legend('topright', c('WT', 'KO'), col=c(red, blue), bty='n',  lwd=1.5, ncol = 2)
dev.off()

