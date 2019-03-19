library(RColorBrewer)
cols <- brewer.pal(n = 3, 'Set1')

cal_zscore <- function(dumps, N, K) {
  kmer <- read.table(paste0('peaks/crosslink.window.', K, '.list'), stringsAsFactors = F)[,1]
  x <- matrix(0, nrow=length(kmer), ncol=101)
  rownames(x) <- kmer
  colnames(x) <- 0:100
  
  input <- read.table(dumps, stringsAsFactors=F)
  x[input[, 1], 1] <- input[, 2]
  
  files <- list.files(paste0('peaks/', K, '/random', N), '*.dump', full.names = T)
  for (i in 1:length(files)) {
    input <- read.table(files[i], stringsAsFactors=F)
    x[input[, 1], i+1] <- input[, 2]
  }
  
  z_score <- (x[, 1]- apply(x[, -1], 1, mean)) / apply(x[, -1], 1, sd)
  names(z_score) <- gsub('T', 'U', names(z_score))
  z_score <- data.frame(z_score=z_score)
  return(z_score)
}

K <- '6mer'

for (N in 1:10) {
  
  zscore1 <- cal_zscore(paste0('peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_', K, '.dumps'), N, K)
  zscore2 <- cal_zscore(paste0('peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_', K, '.dumps'), N, K)
  zscore3 <- cal_zscore(paste0('peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_', K, '.dumps'), N, K)
  zscore4 <- cal_zscore(paste0('peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_', K, '.dumps'), N, K)
  zscore5 <- cal_zscore(paste0('peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_', K, '.dumps'), N, K)
  zscore6 <- cal_zscore(paste0('peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_', K, '.dumps'), N, K)
  zscore <- cbind(zscore1, zscore2, zscore3, zscore4, zscore5, zscore6)
  colnames(zscore) <- paste(colnames(zscore), c('WT1', 'WT2', 'WT3', 'KO1', 'KO2', 'KO3'), sep='.')
  zscore$WT_mean <- (zscore$z_score.WT1 + zscore$z_score.WT2 + zscore$z_score.WT3) / 3
  zscore$KO_mean <- (zscore$z_score.KO1 + zscore$z_score.KO2 + zscore$z_score.KO3) / 3
  zscore$mean_diff <- zscore$KO_mean - zscore$WT_mean
  zscore$ttest.pvalue <- apply(zscore, 1, function(x){t.test(x[1:3],x[4:6])$p.value})
  zscore <- zscore[order(zscore$ttest.pvalue), ]
  
  if (K == '5mer') {
    motifs <- read.table('motif_hnRNP')
    colnames(motifs) <- c('motif', 'class')
    motifs <- merge(motifs, zscore, by.x=1, by.y=0, all.x = T)
    write.table(zscore, paste0('peaks/', K, '/zscore_all_random', N, '.tsv'), quote=F, sep='\t', col.names=NA)
    write.table(motifs[motifs$ttest.pvalue<0.05,], paste0('peaks/', K, '/zscore_motif_hnRNP_random', N, '.tsv'), quote=F, sep='\t', col.names=NA)
    
    # motifs <- motifs[motifs$class != 'Fox',]
    # x <- sort(motifs$rank_diff)
    # n <- motifs$motif[order(motifs$rank_diff)]
    # y <- as.character(motifs$class[order(motifs$rank_diff)])
    # y[y=='green'] <- cols[3]
    # y[y=='blue'] <- cols[2]
    # 
    # pdf(paste0('figure/', K, '_zscore_motif_hnRNP_random', N, '.pdf'))
    # barplot(x, names.arg = n, las=2, col=y, border = NA, main='Z-score rank difference of hnRNP motifs', ylab='WT_specific - KO_specific', space = 0.5, ylim=c(-320, 80), yaxt='n')
    # axis(2, at=seq(-320, 80, 40), lab=seq(-320, 80, 40))
    # legend('bottomright', c('hnRNP M motifs', 'hnRNP C motifs', 'hnRNP H motifs'), fill = c(cols[2], cols[3], 'gray'), bty='n', border = NA)
    # dev.off()
  }
  else {
    write.table(zscore, paste0('peaks/', K, '/zscore_all_random', N, '.tsv'), quote=F, sep='\t', col.names=NA)
  }
}
