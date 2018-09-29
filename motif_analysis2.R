#! /usr/bin/env Rscript

cal_zscore <- function(dumps) {
  kmer <- read.table('merge/crosslink.window.5mer.list')[,1]
  x <- matrix(0, nrow=length(kmer), ncol=101)
  rownames(x) <- kmer
  colnames(x) <- 0:100
  
  input <- read.table(dumps, stringsAsFactors=F)
  x[input[, 1], 1] <- input[, 2]
  
  files <- list.files('merge/random', '*.dump', full.names = T)
  for (i in 1:length(files)) {
    input <- read.table(files[i], stringsAsFactors=F)
    x[input[, 1], i+1] <- input[, 2]
  }
  
  z_score <- (x[, 1]- apply(x[, -1], 1, mean)) / apply(x[, -1], 1, sd)
  names(z_score) <- gsub('T', 'U', names(z_score))
  z_score <- data.frame(score=sort(z_score, decreasing=T), rank=1:length(z_score))
  return(z_score)
}


KO_zscore <- cal_zscore('merge/KO_crosslink_sites_intron_filt_motif_50K_window_5mer.dumps')
WT_zscore <- cal_zscore('merge/WT_crosslink_sites_intron_filt_motif_50K_window_5mer.dumps')
colnames(KO_zscore) <- c('KO_zscore', 'KO_zscore_rank')
colnames(WT_zscore) <- c('WT_zscore', 'WT_zscore_rank')

all <- merge(WT_zscore, KO_zscore, by.x = 0, by.y = 0)
all <- all[order(all$WT_zscore_rank), ]
all$rank_diff <- all$WT_zscore_rank - all$KO_zscore_rank
colnames(all)[1] <- 'motif'

motifs <- read.table('motif_from_HJY')
colnames(motifs) <- c('motif', 'class')
motifs <- merge(motifs, all, by.x=1, by.y=1, all.x = T)
motifs <- motifs[order(motifs$WT_zscore_rank),]
write.table(all, 'motif2/zscore_all.tsv', quote=F, sep='\t', row.names=F)
write.table(motifs, 'motif2/zscore_motif_from_HJY.tsv', quote=F, sep='\t', row.names=F)
