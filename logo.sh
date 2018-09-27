## logo
#wget https://github.com/xuebingwu/kpLogo/archive/v1.1.tar.gz -O kpLogo_v1.1.tar.gz
#tar zxvf kpLogo_v1.1.tar.gz

GENOME=$HOME/Gmatic7/iCount/mus_musculus.88.fa
sort -k5 -n -r merge/merge_ALL_crosslink_sites_sig.bed|awk '{if($5>=10)print}'|awk '{print $1"\t"$2-15"\t"$3+15"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > logo/crosslink.window.bed
bedtools getfasta -fi $GENOME -bed logo/crosslink.window.bed -name -s -fo logo/crosslink.window.fa
kpLogo-1.1/bin/kpLogo logo/crosslink.window.fa -alphabet ACGU -startPos 16 -o logo/kpLogo

