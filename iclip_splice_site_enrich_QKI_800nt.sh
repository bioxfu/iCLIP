module add bedtools/2.25.0
mkdir -p RNAmap/data

awk '{if($5>0) print}' clusters/Mecp2_expt.clusters.bed > RNAmap/data/Mecp2_expt.clusters.filt.bed

cat SE_3wt_3ko_for_RNA_map.csv|awk '{print $1"\t"$3}'|./bin/get_splice_site_800nt.py > RNAmap/data/splice_sig_bin.bed
bedtools intersect -a RNAmap/data/splice_sig_bin.bed -b RNAmap/data/Mecp2_expt.clusters.filt.bed -s -c > RNAmap/data/splice_sig_bin.count.filt
awk '{if($5=="enhanced" && $6=="+" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_en_f.filt
awk '{if($5=="enhanced" && $6=="-" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_en_r.filt
awk '{if($5=="silenced" && $6=="+" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_si_f.filt
awk '{if($5=="silenced" && $6=="-" && $NF>0)print}' RNAmap/data/splice_sig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_sig_rnamap_si_r.filt


cat SE_data_3wt_3ko.csv|awk '{print $1"\tbg"}'|grep 'chr'|sort|uniq|./bin/get_splice_site_800nt.py > RNAmap/data/splice_nonsig_bin.bed
bedtools intersect -a RNAmap/data/splice_nonsig_bin.bed -b RNAmap/data/Mecp2_expt.clusters.filt.bed -s -c > RNAmap/data/splice_nonsig_bin.count.filt
awk '{if($6=="+" && $NF>0)print}' RNAmap/data/splice_nonsig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_nonsig_rnamap_f.filt
awk '{if($6=="-" && $NF>0)print}' RNAmap/data/splice_nonsig_bin.count.filt|sort -k7,7n -k8,8n|bedtools groupby -g 7,8 -c 4 -o count_distinct|sed 's/\t/./' > RNAmap/data/splice_nonsig_rnamap_r.filt


EN_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'enhanced'|wc -l`
SI_EXON=`cut -f4,5 RNAmap/data/splice_sig_bin.bed|sort|uniq|grep 'silenced'|wc -l`
BG_EXON=`awk '{print $4}' RNAmap/data/splice_nonsig_bin.bed|sort|uniq|wc -l`
wc -l RNAmap/data/splice_*sig_rnamap*

Rscript ./RNAmap_800nt.R RNAmap/data $EN_EXON $SI_EXON $BG_EXON 0.05


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


