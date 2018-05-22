FEATURE=$1
GENE=$2

# combine and cluster (tag>1 as filter)
mkdir -p clusters/anno
mkdir -p clusters/stat
echo -e "file\tdata\ttags\tsites\tclusters" > clusters/stat/tag_site_clusters.counts.tsv

for F in $(find sites/*rep*.sites.bed -printf "%f\n"|sed 's/_rep.*//'|sort|uniq); do
    cat sites/${F}_rep*.sites.bed|sort -k4|groupBy -g 4 -c 5 -o sum -full|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6}'|sortBed > clusters/${F}.sites.bed
    
    cat clusters/${F}.sites.bed|bedtools cluster -d 30 > clusters/${F}.clusters.bed
        
    cat clusters/${F}.clusters.bed|cut -f5,7|groupBy -g 2 -c 1,1 -o sum,count > clusters/${F}.clusters.counts
    
    bin/iCLIP_cluster_cutoff.R clusters/${F}.clusters.counts clusters/anno/${F}.clusters.counts.pdf
    
    join -1 7 -2 1 -t $'\t' <(sort -k7 clusters/${F}.clusters.bed) <(cat clusters/${F}.clusters.counts|awk '{if($2>1)print $1}'|sort)|awk '{OFS="\t";print $2,$3,$4,$5,$6,$7,$1}'|sortBed > clusters/${F}.clusters.filtered.bed

    bedtools intersect -a clusters/${F}.clusters.filtered.bed -b $FEATURE -wa -wb |sort -k4,4 -k11,11|groupBy -g 4,5 -c 11 -o distinct|sed -r 's/,.+//'|sed 's/.\.//'> clusters/anno/${F}.clusters.filtered.feature.tsv

    sort -k3 clusters/anno/${F}.clusters.filtered.feature.tsv|groupBy -g 3 -c 2 -o sum |sort -k2nr > clusters/anno/${F}.clusters.filtered.feature.summary.tsv

    bin/iCLIP_cluster_anno.R clusters/anno/${F}.clusters.filtered.feature.summary.tsv clusters/anno/${F}.clusters.filtered.feature.summary.pdf

    cat clusters/${F}.clusters.counts|awk -v f=${F} '{tag+=$2;site+=$3;cluster+=1}END{print f"\tall\t"tag"\t"site"\t"cluster}' >> clusters/stat/tag_site_clusters.counts.tsv

    cat clusters/${F}.clusters.counts|awk -v f=${F} '{if($2>1){tag+=$2;site+=$3;cluster+=1}}END{print f"\tfilt\t"tag"\t"site"\t"cluster}' >> clusters/stat/tag_site_clusters.counts.tsv
    
    cat clusters/${F}.clusters.filtered.bed|sort -k7,7n -k2,2n -k3,3n|groupBy -g 7 -c 2,3,7,5 -o min,max,count,sum -full|cut -f1,8-11 > clusters/${F}.clusters.filtered.grouped.bed
    
    bedtools intersect -a clusters/${F}.clusters.filtered.grouped.bed -b $GENE -wa -wb |cut -f1-5,14|sed -r 's/;.+//'|sed 's/gene_id "//'|sed 's/"//'|sed -r 's/\..+$//'|sortBed|uniq > clusters/${F}.clusters.filtered.grouped.gene.bed
done   

rm clusters/*sites.bed
