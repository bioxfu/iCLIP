TYPE=$1

grep 'Number_of_input_reads' table/*_aligned_stat|sed -r 's/:\s+?/\t/g'|sed 's/_aligned_stat//'|sed 's/table\///' > table/stat.tsv
grep 'Number_of_aligned_reads' table/*_aligned_stat|sed -r 's/:\s+?/\t/g'|sed 's/_aligned_stat//'|sed 's/table\///' >> table/stat.tsv
grep 'Number_of_unique_aligned_reads' table/*_aligned_stat|sed -r 's/:\s+?/\t/g'|sed 's/_aligned_stat//'|sed 's/table\///' >> table/stat.tsv

if [ $TYPE == 'no_rep' ]; then bin/stats_no_rep.R; else bin/stats.R; fi
