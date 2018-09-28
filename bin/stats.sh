TYPE=$1

grep 'Number of input reads' mapping/*/Log.final.out|sed -r 's/:\s+?/\t/'|sed  's/ |//'|sed 's/ /_/g'|sed 's/\/Log.final.out//'|sed 's/mapping\///' > mapping/stat.tsv
grep 'Uniquely mapped reads number' mapping/*/Log.final.out|sed -r 's/:\s+?/\t/'|sed  's/ |//'|sed 's/ /_/g'|sed 's/\/Log.final.out//'|sed 's/mapping\///' >> mapping/stat.tsv

if [ $TYPE == 'no_rep' ]; then bin/stats_no_rep.R; else bin/stats.R; fi
