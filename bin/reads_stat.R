library(fastqcr)
library(magrittr)

qc <- qc_aggregate('fastqc') %>% qc_stats()
qc <- qc[, c('sample', 'tot.seq')]

num <- read.table('stat/all_num')
colnames(num) <- c('sample', 'tot.seq')
num$sample <- sub('.+/', '', sub('.fa', '', num$sample))

qc$clean <- num[grep('.clean', num$sample), 2]
qc$unique_mapped <- num[grep('.map2genome', num$sample), 2] + num[grep('.map2trans', num$sample), 2]

write.table(qc, 'stat/reads_stat.tsv', sep='\t', row.names = F, quote=F)
