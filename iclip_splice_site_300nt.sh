module add bedtools/2.25.0
mkdir -p RNAmap/data

cat mouse_MeCP2_KO_124_SE_RPKM_2_without_NA_enhance_silence|sed 's/chr//'|./bin/get_splice_site_800nt.py > RNAmap/data/splice_sig_bin.bed
cat mouse_MeCP2_KO_124_SE_FDR0.95|sed 's/chr//'|awk '{print $1"\tbg"}'|sort|uniq|./bin/get_splice_site_800nt.py > RNAmap/data/splice_nonsig_bin.bed
EN_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'enhanced'|wc -l`
SI_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'silenced'|wc -l`
BG_EXON=`awk '{print $4}' RNAmap/data/splice_nonsig_bin.bed|sort|uniq|wc -l`

#WT_1=peaks/Mecp2_mouse_rep1.ACATCG_reads_unique_peaks.bedGraph
#WT_2=peaks/Mecp2_mouse_rep2.GCCTAA_reads_unique_peaks.bedGraph
#WT_3=peaks/Mecp2_mouse_rep3.TGGTCA_reads_unique_peaks.bedGraph
#bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_1 -wa -wb > RNAmap/data/WT1_splice_sig_bin_count
#bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_2 -wa -wb > RNAmap/data/WT2_splice_sig_bin_count
#bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_3 -wa -wb > RNAmap/data/WT3_splice_sig_bin_count
#awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_en_f.filt
#awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_en_r.filt
#awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_si_f.filt
#awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_si_r.filt
#awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_en_f.filt
#awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_en_r.filt
#awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_si_f.filt
#awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_si_r.filt
#awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_en_f.filt
#awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_en_r.filt
#awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_si_f.filt
#awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o mean|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_si_r.filt
#Rscript ./RNAmap_300nt.R RNAmap/data/ $EN_EXON $SI_EXON

awk '{if($5>0) print}' merge/merge_ALL_crosslink_sites.bed > RNAmap/data/merge_ALL_crosslink_sites.bed

bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b RNAmap/data/merge_ALL_crosslink_sites.bed -s -c > RNAmap/data/splice_sig_bin.count.filt
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_si_r.filt

bedtools intersect -a RNAmap/data/splice_nonsig_bin.bed -b RNAmap/data/merge_ALL_crosslink_sites.bed -s -c > RNAmap/data/splice_nonsig_bin.count.filt
awk '{if($6=="+" && $NF>0)print}' RNAmap/data/splice_nonsig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_nonsig_rnamap_f.filt
awk '{if($6=="-" && $NF>0)print}' RNAmap/data/splice_nonsig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_nonsig_rnamap_r.filt

Rscript ./RNAmap_800nt.R RNAmap/data/ $EN_EXON $SI_EXON $BG_EXON 0.05
Rscript ./RNAmap_300nt.R RNAmap/data/ $EN_EXON $SI_EXON $BG_EXON 0.05
