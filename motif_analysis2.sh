## In order to compare the Z-scores of motifs bwtween two samples, 
## we subsample the same number of crosslink sites, and generate 
## one permutation backgroud for two samples.

GENOME=$HOME/Gmatic7/iCount/homo_sapiens.88.fa
GENOMESIZE=$HOME/Gmatic7/iCount/homo_sapiens.88.chrLength
INTRON=$HOME/Gmatic7/iCount/homo_sapiens.88.intron.bed
CPU=8
K=6

# exclude WT sites from KO sites
bedtools window -a merge/merge_ALL_crosslink_sites_sig_clusters_KO.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_WT.bed -w 40 -v > merge/KO_crosslink_sites_exc_WT.bed
bedtools window -a merge/merge_ALL_crosslink_sites_sig_clusters_WT.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_KO.bed -w 40 -v > merge/WT_crosslink_sites_exc_KO.bed

# intronic sites
bedtools intersect -a merge/KO_crosslink_sites_exc_WT.bed -b $INTRON -wa -u > merge/KO_crosslink_sites_exc_WT_intron.bed
bedtools intersect -a merge/WT_crosslink_sites_exc_KO.bed -b $INTRON -wa -u > merge/WT_crosslink_sites_exc_KO_intron.bed

# filter sites
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' merge/KO_crosslink_sites_exc_WT_intron.bed merge/KO_crosslink_sites_exc_WT_intron_filt_motif.bed
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' merge/WT_crosslink_sites_exc_KO_intron.bed merge/WT_crosslink_sites_exc_KO_intron_filt_motif.bed

# keep same sites
shuf merge/KO_crosslink_sites_exc_WT_intron_filt_motif.bed|head -50000 > merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K.bed
shuf merge/WT_crosslink_sites_exc_KO_intron_filt_motif.bed|head -50000 > merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K.bed

cat merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window.bed
cat merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window.bed

bedtools getfasta -fi $GENOME -bed merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window.bed -name -s -fo merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window.fa
bedtools getfasta -fi $GENOME -bed merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window.bed -name -s -fo merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window.fa

jellyfish count -m $K -s 100M -t $CPU merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window.fa -o merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window.fa -o merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window_${K}mer.jf

jellyfish dump merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window_${K}mer.jf_0 -c -t|sort > merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window_${K}mer.dumps
jellyfish dump merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window_${K}mer.jf_0 -c -t|sort > merge/WT_crosslink_sites_exc_KO_intron_filt_motif_50K_window_${K}mer.dumps

cat merge/*${K}mer.dumps|cut -f1|sort|uniq > merge/crosslink.window.${K}mer.list

rm merge/*.jf_0

## generate 10 random background
mkdir merge/${K}mer
for N in {1..10}
do
mkdir merge/${K}mer/random
for i in {1..100}
do
	echo $i
    bedtools shuffle -i merge/KO_crosslink_sites_exc_WT_intron_filt_motif_50K_window.bed -g $GENOMESIZE -incl $INTRON > merge/${K}mer/random/random.$i.bed
    bedtools getfasta -fi $GENOME -bed merge/${K}mer/random/random.$i.bed -name -s -fo merge/${K}mer/random/random.$i.fa
    jellyfish count -m $K -s 100M merge/${K}mer/random/random.$i.fa -o merge/${K}mer/random/random.${K}mer.$i.jf
    jellyfish dump merge/${K}mer/random/random.${K}mer.$i.jf_0 -c -t |sort > merge/${K}mer/random/random.${K}mer.$i.dumps
	rm merge/${K}mer/random/random.$i.bed merge/${K}mer/random/random.$i.fa merge/${K}mer/random/random.${K}mer.$i.jf_0
done
mv merge/${K}mer/random merge/${K}mer/random${N}
done


