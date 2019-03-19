dfm <- read.table('/data1/HJY/Project/20181121JY/RNA-Seq/table_8samples/RPKM_table_FDR0.05_FC1.5_all.tsv', header = T, row.names = 1)[1:8]
all <- read.table('all_SE')
un <- read.table('MeCP2_KO_RBFOX2_KI_SE_RPKM_2_union')$V1

write.table(setdiff(all[all$V2 %in% rownames(dfm)[rowSums(dfm) / 8 > 2], 1], un), 'MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement', row.names = F, col.names = F, quote = F)
