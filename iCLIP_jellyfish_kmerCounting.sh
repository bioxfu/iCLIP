WINDOW=21nt.window
K=6
CPU=8

GENOME=$HOME/Gmatic6/genome/mm10/mm10.fa
GENOMESIZE=$HOME/Gmatic6/genome/mm10/mm10.genome.size
EXONS=$HOME/Gmatic6/gene/mm10_vM17/gencode.vM17.exons.slim.bed
NAME=Mecp2

mkdir -p motif/$WINDOW
#sort -k5 -n -r clusters/${NAME}_expt.clusters.filtered.bed|head -1000000 > motif/$WINDOW/crosslink.bed
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

