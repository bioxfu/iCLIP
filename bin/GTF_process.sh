zcat /home/xfu/Gmatic6/gene/mm10_vM17/gencode.vM17.annotation.gtf.gz|awk '{if($3=="exon"){OFS="\t";print $1,$4-1,$5,$12,$5-$4+1,$7,$14,$22}}'|sed 's/[;"]//g' > temp1
sort -k4 temp1|groupBy -g 4 -c 4 -o count > temp2
join -1 4 -2 1 <(sort -k4 temp1) <(sort -k1 temp2)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$7,$8,$9}'|sort -k4,4 -k8,8n > /home/xfu/Gmatic6/gene/mm10_vM17/gencode.vM17.exons.bed
bin/GTF_trans_exons2genome.R /home/xfu/Gmatic6/gene/mm10_vM17/gencode.vM17.exons.bed /home/xfu/Gmatic6/gene/mm10_vM17/gencode.vM17.exons2genome.bed 
cat /home/xfu/Gmatic6/gene/mm10_vM17/gencode.vM17.exons.bed |awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > /home/xfu/Gmatic6/gene/mm10_vM17/gencode.vM17.exons.slim.bed


NAME=vM4

zcat gencode.$NAME.annotation.gtf.gz|awk '{if($3=="exon"){OFS="\t";print $1,$4-1,$5,$12,$5-$4+1,$7,$14,$26}}'|sed 's/[;"]//g' > temp1
sort -k4 temp1|groupBy -g 4 -c 4 -o count > temp2 
join -1 4 -2 1 <(sort -k4 temp1) <(sort -k1 temp2)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$7,$8,$9}'|sort -k4,4 -k8,8n > gencode.$NAME.exons.bed

GTF_trans_exons2genome.r gencode.$NAME.exons.bed gencode.$NAME.exons2genome.bed

rm temp*

#cat $NAME.exons2genome.bed|awk '{if($8==1)print}'|sort -k1,1 -k2,2n|bedtools merge -nms|cut -f4|sed -r 's/;.+//'|tr '|' '\t'> $OUTPUT.first.exons.bed
#cat $NAME.exons2genome.bed|awk '{if($8==$9)print}'|sort -k1,1 -k2,2n|bedtools merge -nms|cut -f4|sed -r 's/;.+//'|tr '|' '\t'> $OUTPUT.last.exons.bed


#cat $OUTPUT.exons.bed|./exon_junction.py > $OUTPUT.exons.junc.bed
#cat $OUTPUT.exons.bed |awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > $OUTPUT.exons.slim.bed
#cat $OUTPUT.exons.bed|awk '{if($7==1){print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}}'|sort -k1,1 -k2,2n|bedtools merge -nms|cut -f4|sed -r 's/;.+//'|tr '|' '\t'> $OUTPUT.first.exons.bed
#cat $OUTPUT.exons.bed|awk '{if($7==$8){print $1"\t"$2"\t"$3"\t"$1"|"$2"|"$3}}'|sort -k1,1 -k2,2n|bedtools merge -nms|cut -f4|sed -r 's/;.+//'|tr '|' '\t'> $OUTPUT.last.exons.bed
#rm $OUTPUT.temp*
#grep -v '#' $GTF|awk '{print "chr"$0}'|sed 's/^chrMT/chrM/' > ${GTF}2
#awk '{if($3=="gene") print}' $GTF|sed -r 's/.+gene_id "//'|sed -r 's/"; gene_name "/\t/'|sed -r 's/"; gene_source.+gene_biotype "/\t/'|sed 's/";//'|sort|uniq > $OUTPUT.ID2Name
#awk '{if($3=="gene") print}' ${GTF}2 > ${GTF}2.gene

