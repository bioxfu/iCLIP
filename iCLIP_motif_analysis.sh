WINDOW=21nt.window
K=6

cat motif/$WINDOW/*.dumps|cut -f1|sort|uniq > motif/$WINDOW/crosslink.window.${K}mer.list

bin/iCLIP_motif_enrich.R motif/$WINDOW crosslink.window.${K}mer.list crosslink.window.${K}mer.dumps *random.${K}mer.*.dumps motif/$WINDOW/motifs.${K}mer

grep -v 'UUU' motif/$WINDOW/motifs.${K}mer.zscore.tsv > motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv

bin/iCLIP_plot_zscore.R motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv motif/$WINDOW/motifs.${K}mer.zscore.filtUUU.pdf

rm motif/$WINDOW/crosslink.window.${K}mer.list motif/$WINDOW/crosslink.window.random.${K}mer.*.dumps motif/$WINDOW/crosslink.window.${K}mer.dumps motif/$WINDOW/crosslink.window.bed
