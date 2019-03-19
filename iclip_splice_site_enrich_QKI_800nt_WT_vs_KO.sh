module add bedtools/2.25.0
mkdir -p RNAmap/data

WT_BED=merge/merge_ALL_crosslink_sites_sig_WT.bed
KO_BED=merge/merge_ALL_crosslink_sites_sig_KO.bed

cat MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|sed 's/chr//'|./bin/get_splice_site_800nt.py > RNAmap/data/splice_sig_bin.bed
cat MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement|sed 's/:/|/g'|sed 's/chr//'|awk '{print $1"\tbg"}'|sort|uniq|./bin/get_splice_site_800nt.py > RNAmap/data/splice_nonsig_bin.bed
EN_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'enhanced'|wc -l`
SI_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'silenced'|wc -l`
BG_EXON=`awk '{print $4}' RNAmap/data/splice_nonsig_bin.bed|sort|uniq|wc -l`
wc -l RNAmap/data/splice_*sig_rnamap*

bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $WT_BED -s -c > RNAmap/data/WT_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_nonsig_bin.bed -b $WT_BED -s -c > RNAmap/data/WT_splice_nonsig_bin_count
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/WT_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/WT_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/WT_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/WT_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/WT_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/WT_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/WT_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/WT_splice_sig_rnamap_si_r.filt
awk '{if($6=="+" && $NF>0)print}' RNAmap/data/WT_splice_nonsig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/WT_splice_nonsig_rnamap_f.filt
awk '{if($6=="-" && $NF>0)print}' RNAmap/data/WT_splice_nonsig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/WT_splice_nonsig_rnamap_r.filt

bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b $KO_BED -s -c > RNAmap/data/KO_splice_sig_bin_count
bedtools intersect -a RNAmap/data/splice_nonsig_bin.bed -b $KO_BED -s -c > RNAmap/data/KO_splice_nonsig_bin_count
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/KO_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/KO_splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/KO_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/KO_splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/KO_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/KO_splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/KO_splice_sig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/KO_splice_sig_rnamap_si_r.filt
awk '{if($6=="+" && $NF>0)print}' RNAmap/data/KO_splice_nonsig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/KO_splice_nonsig_rnamap_f.filt
awk '{if($6=="-" && $NF>0)print}' RNAmap/data/KO_splice_nonsig_bin_count|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/KO_splice_nonsig_rnamap_r.filt

Rscript ./RNAmap_800nt.R RNAmap/data/WT_ $EN_EXON $SI_EXON $BG_EXON 0.005
Rscript ./RNAmap_800nt.R RNAmap/data/KO_ $EN_EXON $SI_EXON $BG_EXON 0.005


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


