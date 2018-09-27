#! /usr/bin/env Rscript
library(fastqcr)
library(magrittr)
library(yaml)
config <- yaml.load_file('config.yaml')

qc_raw <- qc_aggregate('fastqc/raw') %>% qc_stats()
qc_clean <- qc_aggregate('fastqc/trim') %>% qc_stats()
qc <- cbind(qc_raw, qc_clean)
qc <- qc[, c(1,4,9)]
colnames(qc) <- c('sample', 'total.reads', 'trimed.adapter.reads')

stat_tab <- matrix(nrow=length(config$samples), ncol=6)
rownames(stat_tab) <- config$samples
colnames(stat_tab) <- c('PCR.dedupliated.reads', 'uniquely.mapped.reads', 'crosslink.sites', 'crosslink.sites.all', 'crosslink.sites.FDR<0.05', 'crosslink.sites.clusters')
mapping_stat <- read.table('mapping/stat.tsv', stringsAsFactors = F)

for (i in 1:length(config$samples)) {
  x <- read.table(paste0('xlsites/', config$samples[i], '_reads_unique.bed'))
  stat_tab[i, 1:3] <- c(mapping_stat[mapping_stat$V1 == config$samples[i], 3], nrow(x))
}

sites_all <- as.numeric(sub(' .+', '', system('wc -l merge/merge_ALL_crosslink_sites.bed', intern = T)))
sites_sig <- as.numeric(sub(' .+', '', system('wc -l merge/merge_ALL_crosslink_sites_sig.bed', intern = T)))
sites_cluster <- as.numeric(sub(' .+', '', system('cut -f7 merge/merge_ALL_crosslink_sites_sig_clusters.bed|sort|uniq|wc -l', intern = T)))

stat_tab[4, 4:6] <- c(sites_all, sites_sig, sites_cluster)
write.table(cbind(qc, stat_tab), 'table/stats_tab.tsv', row.names = F, quote = F, sep = '\t')
