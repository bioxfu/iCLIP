INPUT=$1
OUTPUT=${INPUT}Graph

N=`cat $INPUT |awk 'BEGIN{v=0}{v+=$5}END{print v/1000000}'`

cat $INPUT|awk -v n="$N" '{if($6=="-"){$5=$5*-1} if($5!=0)print "chr"$1"\t"$2"\t"$3"\t"$5/n}' > $OUTPUT
