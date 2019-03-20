module add bedtools/2.25.0
mkdir -p RNAmap/data

#WT_BED=merge/merge_ALL_crosslink_sites_sig_WT.bed
#KO_BED=merge/merge_ALL_crosslink_sites_sig_KO.bed

WT_1=peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph
WT_2=peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph
WT_3=peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph
KO_1=peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph
KO_2=peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph
KO_3=peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph

cat MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|sed 's/chr//'|./bin/get_splice_site_800nt.py > RNAmap/data/splice_sig_bin.bed
EN_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'enhanced'|wc -l`
SI_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'silenced'|wc -l`

bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_1 -wa -wb > RNAmap/data/WT1_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_2 -wa -wb > RNAmap/data/WT2_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_3 -wa -wb > RNAmap/data/WT3_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $KO_1 -wa -wb > RNAmap/data/KO1_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $KO_2 -wa -wb > RNAmap/data/KO2_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $KO_3 -wa -wb > RNAmap/data/KO3_splice_sig_bin_count


awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT1_splice_sig_rnamap_si_r.filt
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/KO1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO1_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/KO1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO1_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/KO1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO1_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/KO1_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO1_splice_sig_rnamap_si_r.filt

awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT2_splice_sig_rnamap_si_r.filt
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/KO2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO2_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/KO2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO2_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/KO2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO2_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/KO2_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO2_splice_sig_rnamap_si_r.filt

awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/WT3_splice_sig_rnamap_si_r.filt
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/KO3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO3_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/KO3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO3_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/KO3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO3_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/KO3_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 12 -o sum|sed 's/\t/./' > RNAmap/data/KO3_splice_sig_rnamap_si_r.filt

Rscript ./RNAmap_800nt_v2.R RNAmap/data/ $EN_EXON $SI_EXON


#cat $OUTPUT/TableS1.csv|awk '{print $1"\t"$3}'|./get_intron_size.py > $OUTPUT/splice_sig_intron_size
#cat $AS/WT.KO.fisher.test.nonsig.SE.txt|awk '{print $1"\t"$NF}'|grep 'chr'|sort|uniq|./get_intron_size.py > $OUTPUT/splice_nonsig_intron_size
###
#cat TableS1.csv|awk '{print $1"\t"$3}'|./script/get_intron_300nt.py > $OUTPUT/intron_300nt.bed
#bedtools getfasta -fi /home/xfu/Gmatic5/genome/hg19/hg19.fa -bed $OUTPUT/intron_300nt.bed -tab -s -name -fo $OUTPUT/intron_300nt.tab
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced|upstream'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced|downstream'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced'|sed -r 's/:.+//'|sort|uniq -d|wc -l
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced'|sed -r 's/:.+//'|sort|uniq -d > $OUTPUT/olp_event_silenced
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced|5end'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced|3end'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced'|sed -r 's/:.+//'|sort|uniq -d|wc -l
#cat $OUTPUT/intron_300nt.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced'|sed -r 's/:.+//'|sort|uniq -d > $OUTPUT/olp_event_enhanced
#cat TableS1.csv|awk '{print $1"\t"$3}'|./script/get_intron.py > $OUTPUT/intron.bed
#bedtools getfasta -fi /home/xfu/Gmatic5/genome/hg19/hg19.fa -bed $OUTPUT/intron.bed -tab -s -name -fo $OUTPUT/intron.tab
#cat $OUTPUT/intron.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced|upstream'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced|downstream'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'silenced'|sed -r 's/:.+//'|sort|uniq -d|wc -l
#cat $OUTPUT/intron.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced|upstream'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced|downstream'|sed -r 's/:.+//'|sort|uniq|wc -l
#cat $OUTPUT/intron.tab|grep -i 'ACTAA[CT]'|cut -f1|grep 'enhanced'|sed -r 's/:.+//'|sort|uniq -d|wc -l


