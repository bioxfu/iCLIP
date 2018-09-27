INPUT=$1
OUTPUT=$2

echo -e 'type\tsites\treads' > $OUTPUT

grep    'CDS' $INPUT > $INPUT.CDS
grep -v 'CDS' $INPUT|grep 'UTR5' > $INPUT.UTR5
grep -v 'CDS' $INPUT|grep -v 'UTR5'|grep 'UTR3' > $INPUT.UTR3
grep -v 'CDS' $INPUT|grep -v 'UTR5'|grep -v 'UTR3'|grep 'ncRNA' > $INPUT.ncRNA
grep -v 'CDS' $INPUT|grep -v 'UTR5'|grep -v 'UTR3'|grep -v 'ncRNA'|grep 'intron' > $INPUT.intron
grep -v 'CDS' $INPUT|grep -v 'UTR5'|grep -v 'UTR3'|grep -v 'ncRNA'|grep -v 'intron'|grep 'intergenic' > $INPUT.intergenic

wc -l $INPUT.*|grep -v 'total'|awk '{print $2"\t"$1}'|sed -r 's/.+\.//' > $OUTPUT.tmp

for file in `find $INPUT.*`
do
	cat $file|awk -F '\t' 'BEGIN {x=0} {x+=$5} END {print x}' >> $OUTPUT.tmp2
done

paste $OUTPUT.tmp $OUTPUT.tmp2 >> $OUTPUT

rm $OUTPUT.tmp*
