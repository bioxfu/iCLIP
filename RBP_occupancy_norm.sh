GENE=/home/xfu/Gmatic7/iCount/homo_sapiens.88.gtf.gene.bed
GENE_LEN=/home/xfu/Gmatic7/iCount/homo_sapiens.88.gtf.geneLen
mkdir occupancy

#chr10|+|101232253|101232434|101229610|101229814|101237938|101238237	enhanced

## link the crosslink sites with genes
bedtools intersect -a xlsites/MeCP2_KO_rep1.TGGTCA_reads_unique.bed -b $GENE -wa -wb -s|awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}'|groupBy -g 1,2,3 -c 4 -o collapse -full|grep -v ','|cut -f1-6 > occupancy/MeCP2_KO_rep1.bed
bedtools intersect -a xlsites/MeCP2_KO_rep2.CACTGT_reads_unique.bed -b $GENE -wa -wb -s|awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}'|groupBy -g 1,2,3 -c 4 -o collapse -full|grep -v ','|cut -f1-6 > occupancy/MeCP2_KO_rep2.bed
bedtools intersect -a xlsites/MeCP2_KO_rep3.ATTGGC_reads_unique.bed -b $GENE -wa -wb -s|awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}'|groupBy -g 1,2,3 -c 4 -o collapse -full|grep -v ','|cut -f1-6 > occupancy/MeCP2_KO_rep3.bed

bedtools intersect -a xlsites/MeCP2_WT_rep1.CGTGAT_reads_unique.bed -b $GENE -wa -wb -s|awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}'|groupBy -g 1,2,3 -c 4 -o collapse -full|grep -v ','|cut -f1-6 > occupancy/MeCP2_WT_rep1.bed
bedtools intersect -a xlsites/MeCP2_WT_rep2.ACATCG_reads_unique.bed -b $GENE -wa -wb -s|awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}'|groupBy -g 1,2,3 -c 4 -o collapse -full|grep -v ','|cut -f1-6 > occupancy/MeCP2_WT_rep2.bed
bedtools intersect -a xlsites/MeCP2_WT_rep3.GCCTAA_reads_unique.bed -b $GENE -wa -wb -s|awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}'|groupBy -g 1,2,3 -c 4 -o collapse -full|grep -v ','|cut -f1-6 > occupancy/MeCP2_WT_rep3.bed

## merge the three replicates
cat occupancy/MeCP2_KO_rep*.bed|sortBed|groupBy -g 1,2,3,4,6 -c 5 -o sum |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > occupancy/merge_MeCP2_KO.bed
cat occupancy/MeCP2_WT_rep*.bed|sortBed|groupBy -g 1,2,3,4,6 -c 5 -o sum |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > occupancy/merge_MeCP2_WT.bed

## 5' and 3' splicing site of silined/enhanced exon
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed
grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$6-50"\t"$6+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed
grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$5-100"\t"$5+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed

# 1: at least 1 (>=1)
## normalized the cDNA counts by total counts of each gene
#python bin/RBP_occupancy_norm.py occupancy/merge_MeCP2_KO.bed 1 > occupancy/merge_MeCP2_KO_norm.bed
#python bin/RBP_occupancy_norm.py occupancy/merge_MeCP2_WT.bed 1 > occupancy/merge_MeCP2_WT_norm.bed

python bin/RBP_occupancy_norm_RPKM.py HEK293T_KO_RPKM.tsv occupancy/merge_MeCP2_KO.bed 1 > occupancy/merge_MeCP2_KO_norm.bed
python bin/RBP_occupancy_norm_RPKM.py HEK293T_WT_RPKM.tsv occupancy/merge_MeCP2_WT.bed 1 > occupancy/merge_MeCP2_WT_norm.bed

## link the normalized crosslink sites with splicing sites
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_occu

## calculate the distance between cross link sites and splicing sites
## average the normalized counts at each position
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o mean > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_dist2occu


#grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed
#grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed
#grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed
#grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed
#grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed
#grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'enhanced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed
#grep '+' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$7-100"\t"$7+50"\t.\t.\t"$2}' >  MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed
#grep '-' MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhance_silence|grep 'silenced'|cut -f1|sed 's/chr//'|tr '|' '\t'|sort|uniq|awk '{print $1"\t"$4-50"\t"$4+100"\t.\t.\t"$2}' >> MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed

#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed -b occupancy/merge_MeCP2_KO_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_occu
#bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS.bed -b occupancy/merge_MeCP2_WT_norm.bed -wa -wb -s > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_occu

#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_KO_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_5SS_WT_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_KO_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_5SS_WT_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_KO_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_3SS_WT_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_KO_dist2occu
#cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_occu |awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_3SS_WT_dist2occu


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
