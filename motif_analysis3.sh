## In order to compare the Z-scores of motifs bwtween two samples, 
## we subsample the same number of crosslink sites, and generate 
## one permutation backgroud for two samples.

GENOME=$HOME/Gmatic7/iCount/homo_sapiens.88.fa
GENOMESIZE=$HOME/Gmatic7/iCount/homo_sapiens.88.chrLength
INTRON=$HOME/Gmatic7/iCount/homo_sapiens.88.intron.bed
CPU=8
K=5

# exclude WT sites from KO sites
bedtools window -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_WT.bed -w 40 -v > peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT.bed
bedtools window -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_WT.bed -w 40 -v > peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT.bed
bedtools window -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_WT.bed -w 40 -v > peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT.bed
bedtools window -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_KO.bed -w 40 -v > peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO.bed
bedtools window -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_KO.bed -w 40 -v > peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO.bed
bedtools window -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bed -b merge/merge_ALL_crosslink_sites_sig_clusters_KO.bed -w 40 -v > peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO.bed

# intronic sites
bedtools intersect -a peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT.bed -b $INTRON -wa -u > peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron.bed
bedtools intersect -a peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT.bed -b $INTRON -wa -u > peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron.bed
bedtools intersect -a peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT.bed -b $INTRON -wa -u > peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron.bed
bedtools intersect -a peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO.bed -b $INTRON -wa -u > peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron.bed
bedtools intersect -a peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO.bed -b $INTRON -wa -u > peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron.bed
bedtools intersect -a peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO.bed -b $INTRON -wa -u > peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron.bed

# filter sites
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron.bed peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif.bed
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron.bed peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif.bed
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron.bed peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif.bed
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron.bed peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif.bed
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron.bed peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif.bed
bin/filter_sites_by_motif.sh $GENOME 'GCATG|TGCAT' peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron.bed peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif.bed

# keep same sites
shuf peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif.bed|head -10000 > peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K.bed
shuf peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif.bed|head -10000 > peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K.bed
shuf peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif.bed|head -10000 > peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K.bed
shuf peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif.bed|head -10000 > peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K.bed
shuf peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif.bed|head -10000 > peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K.bed
shuf peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif.bed|head -10000 > peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K.bed

cat peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed
cat peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed
cat peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed
cat peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.bed
cat peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.bed
cat peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.bed

bedtools getfasta -fi $GENOME -bed peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed -name -s -fo peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.fa
bedtools getfasta -fi $GENOME -bed peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed -name -s -fo peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.fa
bedtools getfasta -fi $GENOME -bed peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed -name -s -fo peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.fa
bedtools getfasta -fi $GENOME -bed peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.bed -name -s -fo peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.fa
bedtools getfasta -fi $GENOME -bed peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.bed -name -s -fo peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.fa
bedtools getfasta -fi $GENOME -bed peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.bed -name -s -fo peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.fa

jellyfish count -m $K -s 100M -t $CPU peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.fa -o peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.fa -o peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.fa -o peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.fa -o peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.fa -o peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window.fa -o peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.jf

jellyfish dump peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.jf_0 -c -t|sort > peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.dumps
jellyfish dump peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.jf_0 -c -t|sort > peaks/MeCP2_KO_rep2_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.dumps
jellyfish dump peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.jf_0 -c -t|sort > peaks/MeCP2_KO_rep3_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window_${K}mer.dumps
jellyfish dump peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.jf_0 -c -t|sort > peaks/MeCP2_WT_rep1_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.dumps
jellyfish dump peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.jf_0 -c -t|sort > peaks/MeCP2_WT_rep2_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.dumps
jellyfish dump peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.jf_0 -c -t|sort > peaks/MeCP2_WT_rep3_reads_unique_peaks_exc_KO_intron_filt_motif_10K_window_${K}mer.dumps

cat peaks/*${K}mer.dumps|cut -f1|sort|uniq > peaks/crosslink.window.${K}mer.list

rm peaks/*.jf_0

## generate 10 random background
mkdir peaks/${K}mer
for N in {1..10}
do
mkdir peaks/${K}mer/random
for i in {1..100}
do
	echo $N $i
    bedtools shuffle -i peaks/MeCP2_KO_rep1_reads_unique_peaks_exc_WT_intron_filt_motif_10K_window.bed -g $GENOMESIZE -incl $INTRON > peaks/${K}mer/random/random.$i.bed
    bedtools getfasta -fi $GENOME -bed peaks/${K}mer/random/random.$i.bed -name -s -fo peaks/${K}mer/random/random.$i.fa
    jellyfish count -m $K -s 100M peaks/${K}mer/random/random.$i.fa -o peaks/${K}mer/random/random.${K}mer.$i.jf
    jellyfish dump peaks/${K}mer/random/random.${K}mer.$i.jf_0 -c -t |sort > peaks/${K}mer/random/random.${K}mer.$i.dumps
	rm peaks/${K}mer/random/random.$i.bed peaks/${K}mer/random/random.$i.fa peaks/${K}mer/random/random.${K}mer.$i.jf_0
done
mv peaks/${K}mer/random peaks/${K}mer/random${N}
done


