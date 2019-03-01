mat <- NULL
argv <- commandArgs(T)
input <- argv[1]

introns <- c('100nt', '300nt', 'intron')
flanks <- c(5, 10)

for (i in 1:3) {
  for (j in 1:2) {
    intron <- introns[i]
    flank <- flanks[j]
    WT_1 <- sum(read.table(paste0(input, '_', intron, '_motif_flank', flank, '_WT_rep1'), stringsAsFactors = F)$V2)
    WT_2 <- sum(read.table(paste0(input, '_', intron, '_motif_flank', flank, '_WT_rep2'), stringsAsFactors = F)$V2)
    WT_3 <- sum(read.table(paste0(input, '_', intron, '_motif_flank', flank, '_WT_rep3'), stringsAsFactors = F)$V2)
    
    KO_1 <- sum(read.table(paste0(input, '_', intron, '_motif_flank', flank, '_KO_rep1'), stringsAsFactors = F)$V2)
    KO_2 <- sum(read.table(paste0(input, '_', intron, '_motif_flank', flank, '_KO_rep2'), stringsAsFactors = F)$V2)
    KO_3 <- sum(read.table(paste0(input, '_', intron, '_motif_flank', flank, '_KO_rep3'), stringsAsFactors = F)$V2)
    
    x <- c(WT_1, WT_2, WT_3)
    y <- c(KO_1, KO_2, KO_3)
    
    p <- t.test(x, y, alternative = 'greater')$p.value
    
    mat <- rbind(mat, c(intron, flank, x, mean(x), y, mean(y), p))
  }
}

colnames(mat) <- c('Region', 'flank', 'WT.rep1', 'WT.rep2', 'WT.rep3', 'WT.mean', 'KO.rep1', 'KO.rep2', 'KO.rep3', 'KO.mean', 'T-test.pvalue(one-side)')

write.table(mat, paste0(input, '_motif_flank_cpm_t-test.tsv'), sep='\t', quote = F, row.names = F)
