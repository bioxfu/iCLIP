source activate icount

if [ ! -d fastq ]; then
	mkdir -p fastq fastqc/raw fastqc/trim trim dedup mapping xlsites peaks reproduce merge motif figure table track 
fi
