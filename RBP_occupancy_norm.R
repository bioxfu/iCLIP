library(RColorBrewer)
library(scales)

col_set <- brewer.pal(n = 8, 'Set1')

smooth_occu <- function(x, step=3) {
  y <- x
  x <- c(x[1:step], x)
  for (i in 1:length(y)) {
    y[i] <- mean(x[i:(i+step)])
  }
  return(y)
}

ko_si_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_dist2occu')
wt_si_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_dist2occu')
ko_en_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_dist2occu')
wt_en_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_dist2occu')
ko_si_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_dist2occu')
wt_si_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_dist2occu')
ko_en_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_dist2occu')
wt_en_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_dist2occu')


pdf('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_v2.pdf')
par(mfrow=c(2, 2))
par(mar=c(5,4,4,1))
plot(wt_en_3ss$V1, smooth_occu(wt_en_3ss$V2, step=10), ylim=c(0,1), type='l', bty='l', las=2, lwd=1.5, col=col_set[1], ylab='Normalized cDNA counts', xlab='', main='enhanced SE events', xaxt='n')
lines(ko_en_3ss$V1, smooth_occu(ko_en_3ss$V2, step=10), col=col_set[2], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, '3SS', 50))
legend('topleft', c('WT', 'KO'), col = col_set[1:2], lwd = 3, bty='n')
abline(h=0)
rect(100, -0.02, 200, 0.02, col='gray', border = NA)
x <- smooth_occu(wt_en_3ss$V2, step=10)
y <- smooth_occu(ko_en_3ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p) < 0.01 & fc > log2(1.5))
end <- which(p.adjust(p) < 0.01 & fc > log2(1.5))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

par(mar=c(5,1,4,4))
plot(wt_en_5ss$V1, smooth_occu(wt_en_5ss$V2, step=10), ylim=c(0,1), type='l', bty='n', yaxt='n', las=2, lwd=1.5, col=col_set[1], ylab='', xlab='', main='', xaxt='n')
lines(ko_en_5ss$V1, smooth_occu(ko_en_5ss$V2, step=10), col=col_set[2], lwd=1.5)
axis(1, at=c(-100,0,50,100,150,200), lab=c('', -50, '5SS', 50, 100, ''))
abline(h=0)
rect(-50, -0.02, 50, 0.02, col='gray', border = NA)
x <- smooth_occu(wt_en_5ss$V2, step=10)
y <- smooth_occu(ko_en_5ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p) < 0.01 & fc > log2(1.5))
end <- which(p.adjust(p) < 0.01 & fc > log2(1.5))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

par(mar=c(5,4,4,1))
plot(wt_si_3ss$V1, smooth_occu(wt_si_3ss$V2, step=10), ylim=c(0,1), type='l', bty='l', las=2, lwd=1.5, col=col_set[3], ylab='Normalized cDNA counts', xlab='', main='silenced SE events', xaxt='n')
lines(ko_si_3ss$V1, smooth_occu(ko_si_3ss$V2, step=10), col=col_set[4], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, '3SS', 50))
legend('topleft', c('WT', 'KO'), col = col_set[3:4], lwd = 3, bty='n')
abline(h=0)
rect(100, -0.02, 200, 0.02, col='gray', border = NA)
x <- smooth_occu(wt_si_3ss$V2, step=10)
y <- smooth_occu(ko_si_3ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p) < 0.01 & fc > log2(1.5))
end <- which(p.adjust(p) < 0.01 & fc > log2(1.5))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

par(mar=c(5,1,4,4))
plot(wt_si_5ss$V1, smooth_occu(wt_si_5ss$V2, step=10), ylim=c(0,1), type='l', bty='n', yaxt='n', las=2, lwd=1.5, col=col_set[3], ylab='', xlab='', main='', xaxt='n')
lines(ko_si_5ss$V1, smooth_occu(ko_si_5ss$V2, step=10), col=col_set[4], lwd=1.5)
axis(1, at=c(-100,0,50,100,150,200), lab=c('', -50, '5SS', 50, 100, ''))
abline(h=0)
rect(-50, -0.02, 50, 0.02, col='gray', border = NA)
x <- smooth_occu(wt_si_5ss$V2, step=10)
y <- smooth_occu(ko_si_5ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p) < 0.01 & fc > log2(1.5))
end <- which(p.adjust(p) < 0.01 & fc > log2(1.5))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

dev.off()

