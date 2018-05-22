source activate gmatic
conda env export > doc/environment.yml

if [ ! -d fastqc ]; then
	mkdir -p fastq fastqc clean/fasta clean/trim3adapt clusters mapping/map2genome mapping/unmap2genome mapping/map2trans motif reproduce sites stat track 
fi
