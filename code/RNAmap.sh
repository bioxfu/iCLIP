## crosslink site ##
awk '{if($5>1){for(i=0;i<$5;i++) print}}' CellLineA/Expt.crosslink.all.bed > Expt.crosslink.tag.bed

## RPKM 
BAM1=$HOME/Gmatic2/Project/HJY/DEG/YB1/Expt1/WT.tophat_out/accepted_hits.bam
BAM2=$HOME/Gmatic2/Project/HJY/DEG/YB1/Expt1/KO.tophat_out/accepted_hits.bam
#samtools view -c $BAM 
#50849923
#samtools view -c $BAM2
#55759720
TOTAL_READS1=50.849923
TOTAL_READS2=55.759720
bedtools coverage -abam $BAM1 -b $GENE_ROOT/hg19/UCSC_hg19_refGene_represent_exon.bed -counts > refGene_exon.rnaseq.WT &
bedtools coverage -abam $BAM2 -b $GENE_ROOT/hg19/UCSC_hg19_refGene_represent_exon.bed -counts > refGene_exon.rnaseq.KO &
cat refGene_exon.rnaseq.WT|sort -k4|groupBy -g 4,5 -c 7,8 -o sum,sum|awk '{OFS="\t";print $1,$2,($4/51)/($3/1000)}' > refGene.rnaseq.WT.RPKM
cat refGene_exon.rnaseq.KO|sort -k4|groupBy -g 4,5 -c 7,8 -o sum,sum|awk '{OFS="\t";print $1,$2,($4/56)/($3/1000)}' > refGene.rnaseq.KO.RPKM
paste refGene.rnaseq.WT.RPKM refGene.rnaseq.KO.RPKM |awk '{OFS="\t";if($1==$4){print $1,$2,$3,$6}}' > refGene.rnaseq.RPKM
cat refGene.rnaseq.RPKM|awk '{if($3>10)print}' > refGene.rnaseq.RPKM.10
rm refGene.rnaseq.WT.RPKM refGene.rnaseq.KO.RPKM

bedtools coverage -abam $BAM1 -b $GENE_ROOT/hg19/gencode.v19.lincRNA.exon.bed -counts > lincRNA.exon.rnaseq.WT &
bedtools coverage -abam $BAM2 -b $GENE_ROOT/hg19/gencode.v19.lincRNA.exon.bed -counts > lincRNA.exon.rnaseq.KO &
cat lincRNA.exon.rnaseq.WT|sort -k4|groupBy -g 4,5 -c 7,8 -o sum,sum|awk '{OFS="\t";print $1,$2,($4/51)/($3/1000)}' > lincRNA.rnaseq.WT.RPKM
cat lincRNA.exon.rnaseq.KO|sort -k4|groupBy -g 4,5 -c 7,8 -o sum,sum|awk '{OFS="\t";print $1,$2,($4/56)/($3/1000)}' > lincRNA.rnaseq.KO.RPKM
paste lincRNA.rnaseq.WT.RPKM lincRNA.rnaseq.KO.RPKM |awk '{OFS="\t";if($1==$4){print $1,$2,$3,$6}}' > lincRNA.rnaseq.RPKM
cat lincRNA.rnaseq.RPKM|awk '{if($3>1)print}' > lincRNA.rnaseq.RPKM.1
rm lincRNA.rnaseq.WT.RPKM lincRNA.rnaseq.KO.RPKM

## CDS and UTR 100bins ##
join -1 4 -2 2 <(sort -k4 $GENE_ROOT/hg19/UCSC_hg19_refGene_represent_CDS.bed) <(sort -k2 refGene.rnaseq.RPKM.10)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$8,$9}'|sortBed>refGene.rnaseq.RPKM.10.CDS.bed
cut_bed_into_bins.py refGene.rnaseq.RPKM.10.CDS.bed 100 > refGene.rnaseq.RPKM.10.CDS.100bins.bed

join -1 4 -2 2 <(sort -k4 $GENE_ROOT/hg19/UCSC_hg19_refGene_represent_5UTR.bed) <(sort -k2 refGene.rnaseq.RPKM.10)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$8,$9}'|sortBed>refGene.rnaseq.RPKM.10.5UTR.bed
cut_bed_into_bins.py refGene.rnaseq.RPKM.10.5UTR.bed 50 > refGene.rnaseq.RPKM.10.5UTR.50bins.bed

join -1 4 -2 2 <(sort -k4 $GENE_ROOT/hg19/UCSC_hg19_refGene_represent_3UTR.bed) <(sort -k2 refGene.rnaseq.RPKM.10)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$8,$9}'|sortBed>refGene.rnaseq.RPKM.10.3UTR.bed
cut_bed_into_bins.py refGene.rnaseq.RPKM.10.3UTR.bed 50 > refGene.rnaseq.RPKM.10.3UTR.50bins.bed

join -1 4 -2 1 <(sort -k4 $GENE_ROOT/hg19/gencode.v19.lincRNA.gene.bed) <(sort lincRNA.rnaseq.RPKM.1)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$8,$9}'|sortBed>lincRNA.rnaseq.RPKM.1.bed
cut_bed_into_bins.py lincRNA.rnaseq.RPKM.1.bed 100 > lincRNA.rnaseq.RPKM.1.100bins.bed

## iCLIP tag profile 
bedtools coverage -a Expt.crosslink.tag.bed -b refGene.rnaseq.RPKM.10.CDS.100bins.bed  -counts > refGene.rnaseq.RPKM.10.CDS.100bins.bed.tag
bedtools coverage -a Expt.crosslink.tag.bed -b refGene.rnaseq.RPKM.10.5UTR.50bins.bed -counts > refGene.rnaseq.RPKM.10.5UTR.50bins.bed.tag
bedtools coverage -a Expt.crosslink.tag.bed -b refGene.rnaseq.RPKM.10.3UTR.50bins.bed -counts > refGene.rnaseq.RPKM.10.3UTR.50bins.bed.tag
bedtools coverage -a Expt.crosslink.tag.bed -b lincRNA.rnaseq.RPKM.10.100bins.bed -counts > lincRNA.rnaseq.RPKM.10.100bins.bed.tag
bedtools coverage -a Expt.crosslink.tag.bed -b lincRNA.rnaseq.RPKM.1.100bins.bed  -counts > lincRNA.rnaseq.RPKM.1.100bins.bed.tag

cat refGene.rnaseq.RPKM.10.CDS.100bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$9,$10,$10/$7/($3-$2)}'|sort -k1,1 -k6,6n > refGene.rnaseq.RPKM.10.CDS.100bins.bed.tag.norm
sort -k6,6n refGene.rnaseq.RPKM.10.CDS.100bins.bed.tag.norm|groupBy -g 6 -c 8 -o mean >refGene.rnaseq.RPKM.10.CDS.100bins.bed.tag.norm.summary

cat refGene.rnaseq.RPKM.10.5UTR.50bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$9,$10,$10/$7/($3-$2)}'|sort -k1,1 -k6,6n > refGene.rnaseq.RPKM.10.5UTR.50bins.bed.tag.norm
sort -k6,6n refGene.rnaseq.RPKM.10.5UTR.50bins.bed.tag.norm|groupBy -g 6 -c 8 -o mean >refGene.rnaseq.RPKM.10.5UTR.50bins.bed.tag.norm.summary

cat refGene.rnaseq.RPKM.10.3UTR.50bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$9,$10,$10/$7/($3-$2)}'|sort -k1,1 -k6,6n > refGene.rnaseq.RPKM.10.3UTR.50bins.bed.tag.norm
sort -k6,6n refGene.rnaseq.RPKM.10.3UTR.50bins.bed.tag.norm|groupBy -g 6 -c 8 -o mean >refGene.rnaseq.RPKM.10.3UTR.50bins.bed.tag.norm.summary

cat lincRNA.rnaseq.RPKM.1.100bins.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$9,$10,$10/$7/($3-$2)}'|sort -k1,1 -k6,6n > lincRNA.rnaseq.RPKM.1.100bins.bed.tag.norm

./plot_mRNA_profile.r
./plot_lincRNA_profile.r


######  start codon and stop codon ###
cat UCSC_hg19_refGene_represent_CDS.bed|awk '{OFS="\t";if($6=="+"){print $1,$2-300,$2+301,$4,$5,$6,"start"}if($6=="-"){print $1,$3-301,$3+300,$4,$5,$6,"start"}}' > cds_start_stop_300bp.bed
cat UCSC_hg19_refGene_represent_CDS.bed|awk '{OFS="\t";if($6=="+"){print $1,$3-301,$3+300,$4,$5,$6,"stop"}if($6=="-"){print $1,$2-300,$2+301,$4,$5,$6,"stop"}}' >> cds_start_stop_300bp.bed
cat cds_start_stop_300bp.bed|sort -k4,4 -k7,7 > temp; mv temp cds_start_stop_300bp.bed

join -1 4 -2 2 <(sort -k4 cds_start_stop_300bp.bed) <(sort -k2 refGene.rnaseq.RPKM.10)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$7,$9,$10}'|sortBed > cds_start_stop_300bp_RPKM.10.bed

bedtools coverage -a Expt.crosslink.tag.bed -b cds_start_stop_300bp_RPKM.10.bed -d > cds_start_stop_300bp_RPKM.10.bed.tag

cat cds_start_stop_300bp_RPKM.10.bed.tag|awk '{OFS="\t";print $4,$5,$6,$7,$8,$9,$10,$11/$8}' > cds_start_stop_300bp_RPKM.10.bed.tag.norm

./plot_cds_start_stop_300bp_profile_C1-3.r 

cut -f1 cds_start_peak_gene_C1.rpkm.tsv|grep -v 'gene' |sort|uniq > cds_start_peak_gene_C1.genes
cut -f1 cds_start_peak_gene_C3.rpkm.tsv|grep -v 'gene' |sort|uniq > cds_start_peak_gene_C3.genes

join -1 4 -2 1 <(sort -k4 UCSC_hg19_refGene_represent_CDS.bed) <(sort cds_start_peak_gene_C1.rpkm.tsv)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$8,$9,$10}'|sortBed > cds_start_peak_gene_C1.rpkm.CDS.bed
join -1 4 -2 1 <(sort -k4 UCSC_hg19_refGene_represent_CDS.bed) <(sort cds_start_peak_gene_C3.rpkm.tsv)|awk '{OFS="\t";print $2,$3,$4,$1,$5,$6,$8,$9,$10}'|sortBed > cds_start_peak_gene_C3.rpkm.CDS.bed

cat cds_start_peak_gene_C3.rpkm.CDS.bed|awk '{OFS="\t";if($6=="+"){print $1,$2-6,$2+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}if($6=="-"){print $1,$3-6,$3+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}}' > cds_start_peak_gene_C3.rpkm.12bp.bed
cat cds_start_peak_gene_C3.rpkm.CDS.bed|awk '{OFS="\t";if($6=="-"){print $1,$2-6,$2+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}if($6=="+"){print $1,$3-6,$3+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}}' > cds_stop_peak_gene_C3.rpkm.12bp.bed

cat cds_start_peak_gene_C1.rpkm.CDS.bed|awk '{OFS="\t";if($6=="+"){print $1,$2-6,$2+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}if($6=="-"){print $1,$3-6,$3+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}}' > cds_start_peak_gene_C1.rpkm.12bp.bed
cat cds_start_peak_gene_C1.rpkm.CDS.bed|awk '{OFS="\t";if($6=="-"){print $1,$2-6,$2+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}if($6=="+"){print $1,$3-6,$3+6,$4"|"$5"|"$7"|"$8"|"$9"|"$6,".",$6}}' > cds_stop_peak_gene_C1.rpkm.12bp.bed


bedtools getfasta -fi $GENOME_ROOT/hg19/hg19.fa -bed cds_start_peak_gene_C3.rpkm.12bp.bed -s -name -tab -fo temp.1
bedtools getfasta -fi $GENOME_ROOT/hg19/hg19.fa -bed cds_stop_peak_gene_C3.rpkm.12bp.bed -s -name -tab -fo temp.2
head -1 cds_start_peak_gene_C3.rpkm.tsv |awk '{print $0"\tstrand\tstart_12bp\tstop_12bp"}'> cds_start_peak_gene_C3.rpkm.12bp.tsv
paste temp.1 temp.2|cut -f1,2,4|tr '|' '\t'|sort -k3nr >> cds_start_peak_gene_C3.rpkm.12bp.tsv

bedtools getfasta -fi $GENOME_ROOT/hg19/hg19.fa -bed cds_start_peak_gene_C1.rpkm.12bp.bed -s -name -tab -fo temp.1
bedtools getfasta -fi $GENOME_ROOT/hg19/hg19.fa -bed cds_stop_peak_gene_C1.rpkm.12bp.bed -s -name -tab -fo temp.2
head -1 cds_start_peak_gene_C1.rpkm.tsv |awk '{print $0"\tstrand\tstart_12bp\tstop_12bp"}'> cds_start_peak_gene_C1.rpkm.12bp.tsv
paste temp.1 temp.2|cut -f1,2,4|tr '|' '\t'|sort -k3nr >> cds_start_peak_gene_C1.rpkm.12bp.tsv

cut -f2 temp.1 > temp.3
cut -f2 temp.2 > temp.4

rm temp*



