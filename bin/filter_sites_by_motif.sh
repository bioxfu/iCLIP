MOTIF=$1
SITES=$2
OUTPUT=$3

GENOME=$HOME/Gmatic7/iCount/mus_musculus.88.fa

cat $SITES|awk '{if($2>=50)print $1"\t"$2-40"\t"$3+40"\t"$1":"$2":"$3":"$4":"$5":"$6"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > $SITES.window.bed
bedtools getfasta -fi $GENOME -bed $SITES.window.bed -name -s -tab |egrep $MOTIF|sed -r 's/\(.+//'|tr ':' '\t' > $OUTPUT
rm $SITES.window.bed
