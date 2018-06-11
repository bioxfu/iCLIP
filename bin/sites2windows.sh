NAME=$1

cat sites/${NAME}.*.sites.bed |awk '{if($5>=10)print $1"\t"$2-10"\t"$3+10"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > homer/${NAME}.sites.window.bed
