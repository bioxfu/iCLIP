#! /usr/bin/env Rscript
library(fastqcr)
library(magrittr)

qc_raw <- qc_aggregate('fastqc/raw') %>% qc_stats()
qc_clean <- qc_aggregate('fastqc/trim') %>% qc_stats()
qc <- cbind(qc_raw, qc_clean)
qc <- qc[, c(1,4,9)]
colnames(qc) <- c('sample', 'total.reads', 'trimed.adapter.reads')

stat_tab <- matrix(nrow=length(qc$sample), ncol=6)
rownames(stat_tab) <- qc$sample
colnames(stat_tab) <- c('PCR.dedupliated.reads', 'mapped.reads', 'uniquely.mapped.reads', 'crosslink.sites', 'crosslink.sites.FDR<0.05', 'crosslink.sites.clusters')
mapping_stat <- read.table('table/stat.tsv', stringsAsFactors = F)

for (i in 1:length(qc$sample)) {
  x <- read.table(paste0('xlsites/', qc$sample[i], '_reads_unique.bed'))
  y <- read.table(paste0('peaks/', qc$sample[i], '_reads_unique_peaks.bed'))
  z <- read.table(paste0('peaks/', qc$sample[i], '_reads_unique_peaks_clusters.bed'))
  if (length(qc$sample) > 1) {
    stat_tab[i, ] <- c(mapping_stat[mapping_stat$V1 == qc$sample[i], 3], nrow(x), nrow(y), z[nrow(z),7])
  }
  else {
    stat_tab[i, ] <- c(mapping_stat[, 2], nrow(x), nrow(y), z[nrow(z),7])
  }
}

write.table(cbind(qc, stat_tab), 'table/stats_table_final.tsv', row.names = F, quote = F, sep = '\t')
