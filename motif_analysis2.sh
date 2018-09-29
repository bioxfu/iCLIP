## In order to compare the Z-scores of motifs bwtween two samples, 
## we subsample the same number of crosslink sites, and generate 
## one permutation backgroud for two samples.

GENOME=$HOME/Gmatic7/iCount/mus_musculus.88.fa
GENOMESIZE=$HOME/Gmatic7/iCount/mus_musculus.88.chrLength
INTRON=$HOME/Gmatic7/iCount/mus_musculus.88.intron.bed
CPU=8
K=5

# merge Fox1,2,3
cat peaks/KO*peaks.bed|sortBed|bedtools groupby -g 1,2,3 -c 4,5,6 -o distinct,sum,distinct > merge/KO_crosslink_sites.bed
cat peaks/WT*peaks.bed|sortBed|bedtools groupby -g 1,2,3 -c 4,5,6 -o distinct,sum,distinct > merge/WT_crosslink_sites.bed

# intronic sites
bedtools intersect -a merge/KO_crosslink_sites.bed -b $INTRON -wa -u > merge/KO_crosslink_sites_intron.bed
bedtools intersect -a merge/WT_crosslink_sites.bed -b $INTRON -wa -u > merge/WT_crosslink_sites_intron.bed

# filter sites
bin/filter_sites_by_motif.sh 'GCATG|TGCAT' merge/KO_crosslink_sites_intron.bed merge/KO_crosslink_sites_intron_filt_motif.bed
bin/filter_sites_by_motif.sh 'GCATG|TGCAT' merge/WT_crosslink_sites_intron.bed merge/WT_crosslink_sites_intron_filt_motif.bed

# keep same sites
shuf merge/KO_crosslink_sites_intron_filt_motif.bed|head -50000 > merge/KO_crosslink_sites_intron_filt_motif_50K.bed
shuf merge/WT_crosslink_sites_intron_filt_motif.bed|head -50000 > merge/WT_crosslink_sites_intron_filt_motif_50K.bed

cat merge/KO_crosslink_sites_intron_filt_motif_50K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > merge/KO_crosslink_sites_intron_filt_motif_50K_window.bed
cat merge/WT_crosslink_sites_intron_filt_motif_50K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > merge/WT_crosslink_sites_intron_filt_motif_50K_window.bed

bedtools getfasta -fi $GENOME -bed merge/KO_crosslink_sites_intron_filt_motif_50K_window.bed -name -s -fo merge/KO_crosslink_sites_intron_filt_motif_50K_window.fa
bedtools getfasta -fi $GENOME -bed merge/WT_crosslink_sites_intron_filt_motif_50K_window.bed -name -s -fo merge/WT_crosslink_sites_intron_filt_motif_50K_window.fa

jellyfish count -m $K -s 100M -t $CPU merge/KO_crosslink_sites_intron_filt_motif_50K_window.fa -o merge/KO_crosslink_sites_intron_filt_motif_50K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU merge/WT_crosslink_sites_intron_filt_motif_50K_window.fa -o merge/WT_crosslink_sites_intron_filt_motif_50K_window_${K}mer.jf

jellyfish dump merge/KO_crosslink_sites_intron_filt_motif_50K_window_${K}mer.jf_0 -c -t|sort > merge/KO_crosslink_sites_intron_filt_motif_50K_window_${K}mer.dumps
jellyfish dump merge/WT_crosslink_sites_intron_filt_motif_50K_window_${K}mer.jf_0 -c -t|sort > merge/WT_crosslink_sites_intron_filt_motif_50K_window_${K}mer.dumps

cat merge/*.dumps|cut -f1|sort|uniq > merge/crosslink.window.${K}mer.list

rm merge/*.jf_0

mkdir merge/random

for i in {1..100}
do
	echo $i

    bedtools shuffle -i merge/WT_crosslink_sites_intron_filt_motif_50K_window.bed -g $GENOMESIZE -incl $INTRON > merge/random/random.$i.bed
    
    bedtools getfasta -fi $GENOME -bed merge/random/random.$i.bed -name -s -fo merge/random/random.$i.fa
    
    jellyfish count -m $K -s 100M merge/random/random.$i.fa -o merge/random/random.${K}mer.$i.jf
    
    jellyfish dump merge/random/random.${K}mer.$i.jf_0 -c -t |sort > merge/random/random.${K}mer.$i.dumps
    
	rm merge/random/random.$i.bed merge/random/random.$i.fa merge/random/random.${K}mer.$i.jf_0
done



