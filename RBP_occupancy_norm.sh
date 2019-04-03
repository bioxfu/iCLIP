GENE_LEN=/home/xfu/Gmatic7/iCount/homo_sapiens.88.gtf.geneLen
mkdir occupancy

#chr10|+|101232253|101232434|101229610|101229814|101237938|101238237	enhanced

grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed

grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed

grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed

grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed


python bin/RBP_occupancy_norm.py $GENE_LEN merge/merge_ALL_crosslink_sites_sig_KO.bed > occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed
python bin/RBP_occupancy_norm.py $GENE_LEN merge/merge_ALL_crosslink_sites_sig_WT.bed > occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed

bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_occu

bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_occu

grep -v ENSG00000226314 occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_occu > tmp; mv tmp occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_occu

cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_dist2occu

cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_dist2occu




#shuf MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement|head -300 > MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_rand
#grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_rand|sed 's/chr//'|tr ':' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS.bed
#grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_rand|sed 's/chr//'|tr ':' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS.bed
#grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_rand|sed 's/chr//'|tr ':' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS.bed
#grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_rand|sed 's/chr//'|tr ':' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS.bed
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS_KO_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS_WT_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS_KO_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS.bed -b occupancy/merge_ALL_crosslink_sites_sig_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS_WT_occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS_KO_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_3SS_WT_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS_KO_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_complement_5SS_WT_dist2occu
