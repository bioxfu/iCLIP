library(RColorBrewer)
cols <- brewer.pal(n = 3, 'Set1')

cal_zscore <- function(dumps, N) {
  kmer <- read.table('merge/crosslink.window.5mer.list', stringsAsFactors = F)[,1]
  x <- matrix(0, nrow=length(kmer), ncol=101)
  rownames(x) <- kmer
  colnames(x) <- 0:100
  
  input <- read.table(dumps, stringsAsFactors=F)
  x[input[, 1], 1] <- input[, 2]
  
  files <- list.files(paste0('merge/random', N), '*.dump', full.names = T)
  for (i in 1:length(files)) {
    input <- read.table(files[i], stringsAsFactors=F)
    x[input[, 1], i+1] <- input[, 2]
  }
  
  z_score <- (x[, 1]- apply(x[, -1], 1, mean)) / apply(x[, -1], 1, sd)
  names(z_score) <- gsub('T', 'U', names(z_score))
  z_score <- data.frame(score=sort(z_score, decreasing=T), rank=1:length(z_score))
  return(z_score)
}

for (N in 1:10) {
  
  zscore1 <- cal_zscore('merge/WT_crosslink_sites_exc_KO_intron_filt_motif_25K_window_5mer.dumps', N)
  zscore2 <- cal_zscore('merge/KO_crosslink_sites_exc_WT_intron_filt_motif_25K_window_5mer.dumps', N)
  colnames(zscore1) <- c('WT_specific.zscore', 'WT_specific.rank')
  colnames(zscore2) <- c('KO_specific.zscore', 'KO_specific.rank')
  
  all <- merge(zscore1, zscore2, by.x = 0, by.y = 0)
  colnames(all)[1] <- 'motif'
  all <- all[order(all$WT_specific.rank), ]
  all$rank_diff <- all$WT_specific.rank - all$KO_specific.rank
  
  motifs <- read.table('motif_hnRNP')
  colnames(motifs) <- c('motif', 'class')
  motifs <- merge(motifs, all, by.x=1, by.y=1, all.x = T)
  motifs <- motifs[order(motifs$WT_specific.rank),]
  write.table(all, paste0('merge/zscore_all_random', N, '.tsv'), quote=F, sep='\t', row.names=F)
  write.table(motifs, paste0('merge/zscore_motif_hnRNP_random', N, '.tsv'), quote=F, sep='\t', row.names=F)
  
  motifs <- motifs[-c(1:2),]
  x <- sort(motifs$rank_diff)
  n <- motifs$motif[order(motifs$rank_diff)]
  y <- as.character(motifs$class[order(motifs$rank_diff)])
  y[y=='green'] <- cols[3]
  y[y=='blue'] <- cols[2]
  
  pdf(paste0('figure/zscore_motif_hnRNP_random', N, '.pdf'))
  barplot(x, names.arg = n, las=2, col=y, border = NA, main='Z-score rank difference of hnRNP motifs', ylab='WT_specific - KO_specific', space = 0.5, ylim=c(-140, 40), yaxt='n')
  axis(2, at=seq(-140, 40, 20), lab=seq(-140, 40, 20))
  legend('bottomright', c('hnRNP M motifs', 'hnRNP C motifs', 'hnRNP H motifs'), fill = c(cols[2], cols[3], 'gray'), bty='n', border = NA)
  dev.off()
  
}
