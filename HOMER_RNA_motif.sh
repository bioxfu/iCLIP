source activate gmatic
mkdir homer

cat merge/merge_ALL_crosslink_sites_sig.bed|awk '{print "chr"$1"\t"$2-10"\t"$3+10"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > homer/merge_ALL_crosslink_sites_sig_21nt_window.bed
cat merge/merge_ALL_crosslink_sites_sig.bed|awk '{print "chr"$1"\t"$2-5"\t"$3+5"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > homer/merge_ALL_crosslink_sites_sig_11nt_window.bed

findMotifsGenome.pl homer/merge_ALL_crosslink_sites_sig_21nt_window.bed mm10 homer/merge_ALL_crosslink_sites_sig_21nt_window_MotifOutput_6mer -rna -bg bg/mm10_RefSeq_gene.bed -len 6
findMotifsGenome.pl homer/merge_ALL_crosslink_sites_sig_21nt_window.bed mm10 homer/merge_ALL_crosslink_sites_sig_21nt_window_MotifOutput_5mer -rna -bg bg/mm10_RefSeq_gene.bed -len 5
findMotifsGenome.pl homer/merge_ALL_crosslink_sites_sig_11nt_window.bed mm10 homer/merge_ALL_crosslink_sites_sig_11nt_window_MotifOutput_6mer -rna -bg bg/mm10_RefSeq_gene.bed -len 6
findMotifsGenome.pl homer/merge_ALL_crosslink_sites_sig_11nt_window.bed mm10 homer/merge_ALL_crosslink_sites_sig_11nt_window_MotifOutput_5mer -rna -bg bg/mm10_RefSeq_gene.bed -len 5

