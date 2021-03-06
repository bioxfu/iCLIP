# reproducibility of cross-link nucleotide positions

echo -e "file\trep\twindow\tc0\tc1\tc2" > reproduce/crosslink_sites_reproduce_sig.tsv
for F in $(find peaks/*rep*unique_peaks.bed|sed 's/_rep.*//'|sort|uniq); do
    cat ${F}_rep1*unique_peaks.bed |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}} END {print f"\t1\t0\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
    cat ${F}_rep[23]*unique_peaks.bed|bedtools window -a ${F}_rep1.*unique_peaks.bed -b - -w 5 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t5\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
    cat ${F}_rep[23]*unique_peaks.bed|bedtools window -a ${F}_rep1.*.bed -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t30\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv

    cat ${F}_rep2*unique_peaks.bed |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}} END {print f"\t2\t0\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
    cat ${F}_rep[13]*unique_peaks.bed|bedtools window -a ${F}_rep2.*unique_peaks.bed -b - -w 5 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t5\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
    cat ${F}_rep[13]*unique_peaks.bed|bedtools window -a ${F}_rep2.*unique_peaks.bed -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t30\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv

    cat ${F}_rep3*unique_peaks.bed |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}} END {print f"\t3\t0\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
    cat ${F}_rep[12]*unique_peaks.bed|bedtools window -a ${F}_rep3.*unique_peaks.bed -b - -w 5 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t5\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
    cat ${F}_rep[12]*unique_peaks.bed|bedtools window -a ${F}_rep3.*unique_peaks.bed -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t30\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_sig.tsv
done

bin/reproduced_sites.R reproduce/crosslink_sites_reproduce_sig.tsv reproduce/crosslink_sites_reproduce_sig.pdf


echo -e "file\trep\twindow\tc0\tc1\tc2" > reproduce/crosslink_sites_reproduce_all.tsv
for F in $(find xlsites/*rep*unique.bed|sed 's/_rep.*//'|sort|uniq); do
    cat ${F}_rep1*unique.bed |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}} END {print f"\t1\t0\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[23]*unique.bed|bedtools window -a ${F}_rep1.*unique.bed -b - -w 5 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t5\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[23]*unique.bed|bedtools window -a ${F}_rep1.*unique.bed -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t1\t30\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv

    cat ${F}_rep2*unique.bed |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}} END {print f"\t2\t0\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[13]*unique.bed|bedtools window -a ${F}_rep2.*unique.bed -b - -w 5 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t5\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[13]*unique.bed|bedtools window -a ${F}_rep2.*unique.bed -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t2\t30\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv

    cat ${F}_rep3*unique.bed |awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}} END {print f"\t3\t0\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[12]*unique.bed|bedtools window -a ${F}_rep3.*unique.bed -b - -w 5 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t5\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
    cat ${F}_rep[12]*unique.bed|bedtools window -a ${F}_rep3.*unique.bed -b - -w 30 -u|awk -v f=${F} 'BEGIN {c1=0;c2=0;c3=0} {if($5>0){c0+=1}if($5>1){c1+=1}if($5>2){c2+=1}}END{print f"\t3\t30\t"c0"\t"c1"\t"c2}' >> reproduce/crosslink_sites_reproduce_all.tsv
done

bin/reproduced_sites.R reproduce/crosslink_sites_reproduce_all.tsv reproduce/crosslink_sites_reproduce_all.pdf
