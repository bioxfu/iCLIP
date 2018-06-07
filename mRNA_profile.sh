wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz

zcat refGene.txt |awk '{if($14=="cmpl" && $15=="cmpl")print $3"\t"$7"\t"$8"\t"$13"\t"$2"\t"$4}'|sort -k1,1 -k2,2n > refGene_CDS.bed
zcat refGene.txt |awk '{if($14=="cmpl" && $15=="cmpl" && ($7-$5) > 50 && $4=="+") print $3"\t"$5"\t"$7"\t"$13"\t"$2"\t"$4}'|sort -k1,1 -k2,2n > refGene_5UTR.bed
zcat refGene.txt |awk '{if($14=="cmpl" && $15=="cmpl" && ($6-$8) > 50 && $4=="+") print $3"\t"$8"\t"$6"\t"$13"\t"$2"\t"$4}'|sort -k1,1 -k2,2n > refGene_3UTR.bed
zcat refGene.txt |awk '{if($14=="cmpl" && $15=="cmpl" && ($7-$5) > 50 && $4=="-") print $3"\t"$5"\t"$7"\t"$13"\t"$2"\t"$4}'|sort -k1,1 -k2,2n >> refGene_3UTR.bed
zcat refGene.txt |awk '{if($14=="cmpl" && $15=="cmpl" && ($6-$8) > 50 && $4=="-") print $3"\t"$8"\t"$6"\t"$13"\t"$2"\t"$4}'|sort -k1,1 -k2,2n >> refGene_5UTR.bed

## crosslink site ##
awk '{if($5>1){for(i=0;i<$5;i++) print}}' clusters/Mecp2_expt.clusters.bed > clusters/Mecp2_expt.clusters.tag.bed

./cut_bed_into_bins.py refGene_CDS.bed 100 > refGene_CDS.100bins.bed
./cut_bed_into_bins.py refGene_3UTR.bed 50 > refGene_3UTR.50bins.bed
./cut_bed_into_bins.py refGene_5UTR.bed 50 > refGene_5UTR.50bins.bed

## iCLIP tag profile 
bedtools coverage -b clusters/Mecp2_expt.clusters.tag.bed -a refGene_CDS.100bins.bed -counts > refGene_CDS.100bins.bed.tag
bedtools coverage -b clusters/Mecp2_expt.clusters.tag.bed -a refGene_5UTR.50bins.bed -counts > refGene_5UTR.50bins.bed.tag
bedtools coverage -b clusters/Mecp2_expt.clusters.tag.bed -a refGene_3UTR.50bins.bed -counts > refGene_3UTR.50bins.bed.tag

cat refGene_CDS.100bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$8/($3-$2)}'|sort -k1,1 -k4,4n > refGene_CDS.100bins.bed.tag.norm
cat refGene_5UTR.50bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$8/($3-$2)}'|sort -k1,1 -k4,4n > refGene_5UTR.50bins.bed.tag.norm
cat refGene_3UTR.50bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$8/($3-$2)}'|sort -k1,1 -k4,4n > refGene_3UTR.50bins.bed.tag.norm

sort -k4,4n refGene_CDS.100bins.bed.tag.norm|groupBy -g 4 -c 6 -o mean > refGene_CDS.100bins.bed.tag.norm.summary
sort -k4,4n refGene_5UTR.50bins.bed.tag.norm|groupBy -g 4 -c 6 -o mean > refGene_5UTR.50bins.bed.tag.norm.summary
sort -k4,4n refGene_3UTR.50bins.bed.tag.norm|groupBy -g 4 -c 6 -o mean > refGene_3UTR.50bins.bed.tag.norm.summary

Rscript plot_mRNA_profile.r

