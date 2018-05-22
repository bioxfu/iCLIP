export PATH=$PWD/bin:$PATH
EXON=$HOME/Gmatic3/Gene/gencode/gencode_v19_protein_coding_exon.bed

cd ..

# exon-exon or exon-intron junction
mkdir junction
for F in $(find clusters/*_expt.clusters.counts -printf "%f\n"|sed 's/.clusters.counts//'|sort|uniq); do
    cat sites/bed/${F}_rep*.genome.bed|awk '{if($6=="+")print $1"\t"$2"\t"$3+$7"\t"$4"\t.\t"$6;if($6=="-")print $1"\t"$2-$7"\t"$3"\t"$4"\t.\t"$6}'|sortBed|uniq > junction/${F}_genome_tag.bed
    
    bedtools intersect -a junction/${F}_genome_tag.bed -b $EXON -wa -wb |awk '{if(!($2>=$8 && $3<=$9))print}'|cut -f4|sort|uniq > junction/${F}_exon_intron_junction
    
    cat sites/bed/${F}_rep*.trans2genome.bed|cut -f4|sort|uniq > junction/${F}_exon_exon_junction
    
    cut -f4,5 clusters/${F}.clusters.filtered.bed|sort|uniq > junction/${F}.clusters.filtered
    
    join -j 1 -t $'\t' <(sort junction/${F}.clusters.filtered) <(sort junction/${F}_exon_intron_junction) > junction/${F}_exon_intron_junction_tag_count
    
    join -j 1 -t $'\t' <(sort junction/${F}.clusters.filtered) <(sort junction/${F}_exon_exon_junction)   > junction/${F}_exon_exon_junction_tag_count

    bedtools intersect -a junction/${F}_genome_tag.bed -b $EXON -wa -wb|awk '{if(($6=="+" && $2<$8 && $3>$8 && $3<$9)||($6=="-" && $2<$9 && $2>$8 && $3>$9))print}'|cut -f4|sort|uniq > junction/${F}_exon_intron_junction_3SS
    
    bedtools intersect -a junction/${F}_genome_tag.bed -b $EXON -wa -wb|awk '{if(($6=="+" && $2<$9 && $2>$8 && $3>$9)||($6=="-" && $2<$8 && $3>$8 && $3<$9 ))print}'|cut -f4|sort|uniq > junction/${F}_exon_intron_junction_5SS

    join -j 1 -t $'\t' <(sort junction/${F}.clusters.filtered) <(sort junction/${F}_exon_intron_junction_3SS) > junction/${F}_exon_intron_junction_3SS_tag_count
    
    join -j 1 -t $'\t' <(sort junction/${F}.clusters.filtered) <(sort junction/${F}_exon_intron_junction_5SS) > junction/${F}_exon_intron_junction_5SS_tag_count

    cat junction/${F}_exon_intron_junction_tag_count| awk '{ SUM += $2} END { print SUM }'
    
    cat junction/${F}_exon_exon_junction_tag_count| awk '{ SUM += $2} END { print SUM }'

    cat clusters/${F}.clusters.filtered.bed| awk '{ SUM += $5} END { print SUM }'
    
    cat junction/${F}_exon_intron_junction_3SS_tag_count| awk '{ SUM += $2} END { print SUM }'
    
    cat junction/${F}_exon_intron_junction_5SS_tag_count| awk '{ SUM += $2} END { print SUM }'
done


