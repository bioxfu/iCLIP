OUTPUT=/home/xfu/Gmatic2/Project/HJY/iCLIP/YB1/mm0/CellLineB/Expt
SAMPLES=(s1 s2 s3)
RANGE=(0 1 2)

ADAPT='TGAGATCGGAAG'
INDEX=$HOME/Gmatic2/Data/base/hsa/hg19/hg19
TRANS=$HOME/Gmatic2/Data/base/hsa/hg19/hg19_trans
EXONS=$HOME/Gmatic2/Data/base/hsa/hg19/gtf/hg19.exons2genome.bed
ANNO=$HOME/Gmatic2/Data/base/hsa/hg19/gtf/hg19.ID2Name
FEATURE=$HOME/Gmatic2/Data/base/hsa/gencode/defining_genomic_regions/gencode_v19_all.bed
GENES=$HOME/Gmatic2/Data/base/hsa/gencode/defining_genomic_regions/gencode.v19.annotation.gene.gtf
CPUNUM=25
IGV=hg19
BIN=$HOME/Gmatic2/Workflow/iCLIP

mkdir $OUTPUT.log
#4.mapping to genome and transcriptome
for i in ${RANGE[*]}
do
    BARCODE=${BARCODES[$i]}
    SAMPLE=${SAMPLES[$i]}
    echo "ID:$SAMPLE($BARCODE)" >> $OUTPUT.log/bowtie.log
    bowtie -p$CPUNUM -f -v0 -k2 -m1 $INDEX $OUTPUT.$SAMPLE.fa $OUTPUT.$SAMPLE.map 2>> $OUTPUT.log/bowtie.log
    bowtie -p$CPUNUM -f -v0 $INDEX $OUTPUT.$SAMPLE.fa $OUTPUT.$SAMPLE.map.2 --un $OUTPUT.$SAMPLE.un --quiet
    bowtie -p$CPUNUM -f -v0 $TRANS $OUTPUT.$SAMPLE.un $OUTPUT.$SAMPLE.map.trans 2>> $OUTPUT.log/bowtie.log
    rm $OUTPUT.$SAMPLE.map.2 $OUTPUT.$SAMPLE.un
done

#5.finding the crosslink positions
for i in ${RANGE[*]}
do
    BARCODE=${BARCODES[$i]}
    SAMPLE=${SAMPLES[$i]}
    cat $OUTPUT.$SAMPLE.map|$BIN/get_crosslink_pos.py > $OUTPUT.$SAMPLE.bed.temp1
    cat $OUTPUT.$SAMPLE.map.trans|awk '{if($2=="+"&&$4>0){print $3"\t"$4-1"\t"$4}}'|sort -k1,1 -k2,2n|uniq -c|awk '{print $2"\t"$3"\t"$4"\t"$1}' > $OUTPUT.$SAMPLE.bed.temp2
    bedtools intersect -a $EXONS -b $OUTPUT.$SAMPLE.bed.temp2 -wa -wb|$BIN/trans_coor_convert.py >> $OUTPUT.$SAMPLE.bed.temp1
    cat $OUTPUT.$SAMPLE.bed.temp1|sort -k1,1 -k2,2n|bedtools groupby -g 1,2 -c 5 -o sum -full|awk '{print $1"\t"$2"\t"$3"\t"$1":"$3":"$6"\t"$7"\t"$6}' > $OUTPUT.$SAMPLE.crosslink.bed
    rm $OUTPUT.$SAMPLE.bed.temp*
done

#7. combine and cluster
cat $OUTPUT.${SAMPLES[0]}.crosslink.bed $OUTPUT.${SAMPLES[1]}.crosslink.bed $OUTPUT.${SAMPLES[2]}.crosslink.bed|sort -k4|bedtools groupby -g 4 -c 5 -o sum -full|sort -k1,1 -k2,2n|bedtools cluster -s -d 31 |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6"\t"$8}'|sort -k1,1 -k2,2n > $OUTPUT.crosslink.all.bed
cat $OUTPUT.crosslink.all.bed|$BIN/bed2wig.py > $OUTPUT.crosslink.all.wig
igvtools toTDF $OUTPUT.crosslink.all.wig $OUTPUT.crosslink.all.tdf $IGV



