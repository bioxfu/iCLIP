WINDOW=21nt.window
K=5
CPU=8

GENOME=$HOME/Gmatic6/genome/mm10/mm10.fa
GENOMESIZE=$HOME/Gmatic6/genome/mm10/mm10.genome.size
EXONS=$HOME/Gmatic6/gene/mm10_vM17/gencode.vM17.exons.slim.bed
NAME=Mecp2
GENE=$HOME/Gmatic6/gene/mm10_vM17/gencode.vM17.annotation.gtf.gz

mkdir -p motif/$WINDOW
#cat sites/${NAME}_expt_rep[123].*.sites.bed | sort -k5 -n -r > motif/$WINDOW/crosslink.bed
#cat sites/${NAME}_expt_rep[13].*.sites.bed | sort -k5 -n -r > motif/$WINDOW/crosslink.bed
sort -k5 -n -r clusters/${NAME}_expt.clusters.filtered.bed > motif/$WINDOW/crosslink.bed

if [ "$WINDOW" == "21nt.window" ]
then
    cat motif/$WINDOW/crosslink.bed|awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > motif/$WINDOW/crosslink.window.bed
elif [ "$WINDOW" == "11nt.window" ]
then
    cat motif.$WINDOW/crosslink.bed|awk '{print $1"\t"$2-5"\t"$3+5"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > motif.$WINDOW/crosslink.window.bed
fi

bedtools getfasta -fi $GENOME -bed motif/$WINDOW/crosslink.window.bed -name -s -fo motif/$WINDOW/crosslink.window.fa

jellyfish count -m $K -s 100M -t $CPU motif/$WINDOW/crosslink.window.fa -o motif/$WINDOW/crosslink.window.${K}mer.jf

jellyfish dump motif/$WINDOW/crosslink.window.${K}mer.jf_0 -c -t|sort > motif/$WINDOW/crosslink.window.${K}mer.dumps

rm motif/$WINDOW/crosslink.window.${K}mer.jf_0 motif/$WINDOW/crosslink.bed

for i in {0..24}
do
    bin/iCLIP_jellyfish_shuffle_paralle.sh $WINDOW $GENOME $GENOMESIZE $EXONS $K $i &     
done

############
cat motif/$WINDOW/*.dumps|cut -f1|sort|uniq > motif/$WINDOW/crosslink.window.${K}mer.list

bin/iCLIP_motif_enrich.R motif/$WINDOW crosslink.window.${K}mer.list crosslink.window.${K}mer.dumps *random.${K}mer.*.dumps motif/$WINDOW/motifs.${K}mer

grep -v 'UUU' motif/$WINDOW/motifs.${K}mer.zscore.tsv > motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv

bin/iCLIP_plot_zscore.R motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.pdf

rm motif/$WINDOW/crosslink.window.${K}mer.list motif/$WINDOW/crosslink.window.random.${K}mer.*.dumps motif/$WINDOW/crosslink.window.${K}mer.dumps motif/$WINDOW/crosslink.window.bed


#####
bedtools intersect -a motif/$WINDOW/crosslink.bed -b $GENE -wa -wb |awk '{if($10=="gene")print $21"\t"$4}'|sed 's/[";]//g'|sort|uniq > motif/$WINDOW/crosslink.gene_list
bedtools getfasta -fi $GENOME -bed motif/$WINDOW/crosslink.window.bed -name -tab -s -fo motif/$WINDOW/crosslink.window.tab
