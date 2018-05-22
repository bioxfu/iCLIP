GENE=$HOME/Gmatic3/Gene/gencode/gencode.v19.annotation.protein_coding.gtf
MIRNA=$HOME/Gmatic3/Gene/hg19/miRBase_v20_hsa.gff3

## miRNA and protein coding gene
mkdir gene 
for F in $(find clusters/*_expt.clusters.counts -printf "%f\n"|sed 's/.clusters.counts//'|sort|uniq); do
    echo -e 'tag\tcounts\tmiRNA_chrom\ttype\tmiRNA_start\tmiRNA_end\tmiRNA_strand\tmiRNA_anno' > gene/${F}.tag_on_miRNA.tsv
    bedtools intersect -a clusters/${F}.clusters.filtered.bed -b $MIRNA -wa -wb -s |awk -F '\t' '{OFS="\t";print $4,$5,$8,$10,$11,$12,$14,$16}' >> gene/${F}.tag_on_miRNA.tsv
    cut -f2,4,8 gene/${F}.tag_on_miRNA.tsv|sort -k3|grep -v 'miRNA_anno'|groupBy -g 2,3 -c 1 -o sum > gene/${F}.tag_on_miRNA_sum.tsv

    bedtools intersect -a clusters/${F}.clusters.filtered.bed -b $GENE -wa -wb -s |cut -f5,16|sed -r 's/;.+//'|sed 's/gene_id "//'|sed 's/"//'|sort -k2|groupBy -g 2 -c 1 -o sum > gene/${F}.tag_on_protein_coding_gene.tsv
done


