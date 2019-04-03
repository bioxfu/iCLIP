library(RColorBrewer)
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

ko_si_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_dist2occu')
wt_si_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_dist2occu')
ko_en_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_dist2occu')
wt_en_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_dist2occu')

pdf('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_downstream.pdf')
par(mfrow=c(2, 2))
plot(wt_en_5ss$V1, smooth_occu(wt_en_5ss$V2, step=2), type='l', lwd=1.5, col=col_set[1], ylab='Normalized crosslink events', xlab='Position relative to 5SS', main='enhanced SE events', xaxt='n')
lines(ko_en_5ss$V1, smooth_occu(ko_en_5ss$V2, step=2), col=col_set[2], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-50, 0, 50, 100))
legend('topleft', c('WT', 'KO'), col = col_set[1:2], lwd = 3, bty='n')

plot(wt_en_3ss$V1, smooth_occu(wt_en_3ss$V2, step=2), type='l', lwd=1.5, col=col_set[1], ylab='Normalized crosslink events', xlab='Position relative to 3SS', main='enhanced SE events', xaxt='n')
lines(ko_en_3ss$V1, smooth_occu(ko_en_3ss$V2, step=2), col=col_set[2], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, 0, 50))
legend('topleft', c('WT', 'KO'), col = col_set[1:2], lwd = 3, bty='n')

plot(wt_si_5ss$V1, smooth_occu(wt_si_5ss$V2, step=2), type='l', lwd=1.5, col=col_set[3], ylab='Normalized crosslink events', xlab='Position relative to 5SS', main='silenced SE events', xaxt='n')
lines(ko_si_5ss$V1, smooth_occu(ko_si_5ss$V2, step=2), col=col_set[4], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-50, 0, 50, 100))
legend('topleft', c('WT', 'KO'), col = col_set[3:4], lwd = 3, bty='n')

plot(wt_si_3ss$V1, smooth_occu(wt_si_3ss$V2, step=2), type='l', lwd=1.5, col=col_set[3], ylab='Normalized crosslink events', xlab='Position relative to 3SS', main='silenced SE events', xaxt='n')
lines(ko_si_3ss$V1, smooth_occu(ko_si_3ss$V2, step=2), col=col_set[4], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, 0, 50))
legend('topleft', c('WT', 'KO'), col = col_set[3:4], lwd = 3, bty='n')
dev.off()


ko_si_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_dist2occu')
wt_si_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_dist2occu')
ko_en_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_dist2occu')
wt_en_5ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_dist2occu')

ko_si_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_dist2occu')
wt_si_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_dist2occu')
ko_en_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_dist2occu')
wt_en_3ss <- read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_dist2occu')

pdf('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_upstream.pdf')
par(mfrow=c(2, 2))
plot(wt_en_5ss$V1, smooth_occu(wt_en_5ss$V2, step=1), type='l', lwd=1.5, col=col_set[1], ylab='Normalized crosslink events', xlab='Position relative to 5SS', main='enhanced SE events', xaxt='n')
lines(ko_en_5ss$V1, smooth_occu(ko_en_5ss$V2, step=1), col=col_set[2], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-50, 0, 50, 100))
legend('topleft', c('WT', 'KO'), col = col_set[1:2], lwd = 3, bty='n')

plot(ko_en_3ss$V1, smooth_occu(ko_en_3ss$V2, step=1), type='l', lwd=1.5, col=col_set[2], ylab='Normalized crosslink events', xlab='Position relative to 3SS', main='enhanced SE events', xaxt='n')
lines(wt_en_3ss$V1, smooth_occu(wt_en_3ss$V2, step=1), col=col_set[1], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, 0, 50))
legend('topleft', c('WT', 'KO'), col = col_set[1:2], lwd = 3, bty='n')

plot(ko_si_5ss$V1, smooth_occu(ko_si_5ss$V2, step=1), type='l', lwd=1.5, col=col_set[4], ylab='Normalized crosslink events', xlab='Position relative to 5SS', main='silenced SE events', xaxt='n')
lines(wt_si_5ss$V1, smooth_occu(wt_si_5ss$V2, step=1), col=col_set[3], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-50, 0, 50, 100))
legend('topleft', c('WT', 'KO'), col = col_set[3:4], lwd = 3, bty='n')

plot(ko_si_3ss$V1, smooth_occu(ko_si_3ss$V2, step=1), type='l', lwd=1.5, col=col_set[4], ylab='Normalized crosslink events', xlab='Position relative to 3SS', main='silenced SE events', xaxt='n')
lines(wt_si_3ss$V1, smooth_occu(wt_si_3ss$V2, step=1), col=col_set[3], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, 0, 50))
legend('topleft', c('WT', 'KO'), col = col_set[3:4], lwd = 3, bty='n')
dev.off()
