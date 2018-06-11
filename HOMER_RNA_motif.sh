mkdir homer

find sites/*.bed|grep -v 'ctrl'|./bin/rush -k "bin/sites2windows.sh  {%@([^.]+)\.}"

find homer/* | ./bin/rush -k "findMotifsGenome.pl {} mm10 {@(.+).sites.*}_MotifOutput -rna -bg bg/mm10_RefSeq_gene.bed -len 6"
