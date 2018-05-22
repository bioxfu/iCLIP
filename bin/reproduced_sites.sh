# reproducibility of cross-link nucleotide positions
mkdir reproduce

echo -e "file\trep\twindow\tc0\tc1\tc2" > reproduce/reproduce.counts.tsv
for F in $(find sites/*rep*.sites.bed|sed 's/_rep.*//'|sort|uniq); do
    cat ${F}_rep1* |awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t0\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
    cat ${F}_rep[23]*|bedtools window -a ${F}_rep1.*.sites.bed -b - -w 5 -u|awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t5\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
    cat ${F}_rep[23]*|bedtools window -a ${F}_rep1.*.sites.bed -b - -w 30 -u|awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t30\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv

    cat ${F}_rep2* |awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t0\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
    cat ${F}_rep[13]*|bedtools window -a ${F}_rep2.*.sites.bed -b - -w 5 -u|awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t5\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
    cat ${F}_rep[13]*|bedtools window -a ${F}_rep2.*.sites.bed -b - -w 30 -u|awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t30\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv

    cat ${F}_rep3* |awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t0\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
    cat ${F}_rep[12]*|bedtools window -a ${F}_rep3.*.sites.bed -b - -w 5 -u|awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t5\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
    cat ${F}_rep[12]*|bedtools window -a ${F}_rep3.*.sites.bed -b - -w 30 -u|awk -v f=${F} '{if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t30\t"c0"\t"c1"\t"c2}' >> reproduce/reproduce.counts.tsv
done

bin/iCLIP_reprod.R reproduce/reproduce.counts.tsv reproduce/reproduce.counts.pdf
