export PATH=$PWD/bin:$PATH

cd ..

## motifs in reads ##
mkdir -p motif/reads

cat sites/bed/CellLineA_expt_rep*genome.bed |grep '+'|cut -f4 > motif/reads/tmp.id
cat sites/bed/CellLineA_expt_rep*genome.bed |grep '-'|cut -f4 >> motif/reads/tmp.id
cat sites/bed/CellLineA_expt_rep*genome.bed |grep '+'|cut -f5 > motif/reads/tmp.seq
cat sites/bed/CellLineA_expt_rep*genome.bed |grep '-'|cut -f5|tr 'ATCG' 'TAGC'|rev >> motif/reads/tmp.seq
paste motif/reads/tmp.id motif/reads/tmp.seq|sort|groupBy -g 1,2 -c 2 -o count > motif/reads/CellLineA_expt_aln_reads
rm motif/reads/tmp.*

cut -f4 clusters/CellLineA_expt.clusters.filtered.bed|sort|uniq > motif/reads/CellLineA_expt.clusters.filtered
join -j 1 -t $'\t' <(sort motif/reads/CellLineA_expt.clusters.filtered) <(sort motif/reads/CellLineA_expt_aln_reads) > motif/reads/CellLineA_expt_aln_reads_clusters.filtered

head -10 motif/21nt.window/motifs.6mer.zscore.tsv |cut -f1|tr 'U' 'T' > motif/21nt.window/motifs.6mer.top10
head -10 motif/21nt.window/motifs.6mer.zscore.filtUUU.tsv |cut -f1|tr 'U' 'T' > motif/21nt.window/motifs.6mer.filtUUU.top10

iCLIP_motif_search.py reads/CellLineA_expt_aln_reads_clusters.filtered motif/21nt.window/motifs.6mer.filtUUU.top10 
iCLIP_motif_search.py reads/CellLineA_expt_aln_reads_clusters.filtered motif/21nt.window/motifs.6mer.top10 

iCLIP_motif_search_CAUC.py reads/CellLineA_expt_aln_reads_clusters.filtered 


