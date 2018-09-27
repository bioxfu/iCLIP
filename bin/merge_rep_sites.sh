#cat peaks/*expt*.bed|sortBed|bedtools groupby -g 1,2,3 -c 4,5,6 -o distinct,sum,distinct > merge/merge_sig_crosslink_sites.bed
cat xlsites/*expt*unique.bed|sortBed|bedtools groupby -g 1,2,3 -c 4,5,6 -o distinct,sum,distinct > merge/merge_ALL_crosslink_sites.bed
