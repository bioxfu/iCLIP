K=6
WINDOW=11nt.window

GENOME=$HOME/Gmatic7/iCount/mus_musculus.88.fa
GENOMESIZE=$HOME/Gmatic7/iCount/mus_musculus.88.chrLength
EXONS=$HOME/Gmatic7/iCount/mus_musculus.88.pc.gene.bed
CPU=8

mkdir -p motif/$WINDOW
sort -k5 -n -r merge/merge_ALL_crosslink_sites_sig.bed|awk '{if($5>0)print}' > motif/$WINDOW/crosslink.bed

if [ "$WINDOW" == "21nt.window" ]
then
    cat motif/$WINDOW/crosslink.bed|awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > motif/$WINDOW/crosslink.window.bed
elif [ "$WINDOW" == "11nt.window" ]
then
    cat motif/$WINDOW/crosslink.bed|awk '{print $1"\t"$2-5"\t"$3+5"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > motif/$WINDOW/crosslink.window.bed
fi

bedtools getfasta -fi $GENOME -bed motif/$WINDOW/crosslink.window.bed -name -s -fo motif/$WINDOW/crosslink.window.fa
jellyfish count -m $K -s 100M -t $CPU motif/$WINDOW/crosslink.window.fa -o motif/$WINDOW/crosslink.window.${K}mer.jf
jellyfish dump motif/$WINDOW/crosslink.window.${K}mer.jf -c -t|sort > motif/$WINDOW/crosslink.window.${K}mer.dumps

# permutation 100 times
for i in {0..24}; do  echo $i ; done|parallel --gnu "bin/iCLIP_jellyfish_shuffle_paralle.sh $WINDOW $GENOME $GENOMESIZE $EXONS $K {}"

cat motif/$WINDOW/*.dumps|cut -f1|sort|uniq > motif/$WINDOW/crosslink.window.${K}mer.list
bin/iCLIP_motif_enrich.R motif/$WINDOW crosslink.window.${K}mer.list crosslink.window.${K}mer.dumps *random.${K}mer.*.dumps motif/$WINDOW/motifs.${K}mer
grep -v 'UUU' motif/$WINDOW/motifs.${K}mer.zscore.tsv > motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv
bin/iCLIP_plot_zscore.R motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.pdf

rm motif/$WINDOW/crosslink.window.${K}mer.list motif/$WINDOW/crosslink.window.random.${K}mer.*.dumps motif/$WINDOW/crosslink.window.${K}mer.dumps motif/$WINDOW/crosslink.window.bed


