configfile: "config.yaml"

rule all:
	input:
		expand('fastq/{sample}.fastq.gz', sample=config['samples']),
		expand('fastqc/{sample}_fastqc.html', sample=config['samples']),
		expand('clean/{sample}.clean.fa', sample=config['samples']),
		'stat/reads_stat.tsv',
		expand('sites/{sample}.sites.bed', sample=config['samples']),
		expand('track/{sample}.tdf', sample=config['samples']),
		'reproduce/reproduce.counts.tsv',
		'reproduce/reproduce.counts.pdf',
		'clusters/stat/tag_site_clusters.counts.tsv',


rule split_fastq:
	input:
		config['raw_reads']
	output:
		'fastq/{sample}.fastq.gz'
	shell:
		'bin/split_fastq_by_barcode.py {input} {output}'

rule fastqc:
	input:
		'fastq/{sample}.fastq.gz'
	output:
		'fastqc/{sample}_fastqc.html'
	shell:
		'fastqc -o fastqc {input}'

rule fastq2fasta:
	input:
		'fastq/{sample}.fastq.gz'
	output:
		'clean/fasta/{sample}.fa.gz'
	shell:
		"zcat {input} | fastq_to_fasta -Q33 -r -n -v -z -o {output}"

rule trim3adapter:
	input:
		'clean/fasta/{sample}.fa.gz'
	output:
		'clean/trim3adapt/{sample}.trim.fa.gz'
	params:
		adapt = config['adapt']
	shell:
		 "cutadapt -f fasta -a {params.adapt} -m15 --trimmed-only {input} | gzip -c > {output}"

rule remove_PCR_duplicates:
	input:
		'clean/trim3adapt/{sample}.trim.fa.gz'
	output:
		'clean/{sample}.clean.fa'
	shell:
		 "zcat {input} | grep -v '>' | bin/split_barcode.py | sort | uniq | cut -f1 | bin/seq2fasta.py > {output}"

rule reads_map2genome:
	input:
		'clean/{sample}.clean.fa'
	output:
		map2genome = 'mapping/map2genome/{sample}.map2genome'
	params:
		DNA = config['DNA']
	shell:
		 "bowtie -p5 -f -v1 -k2 -m1 {params.DNA} {input} {output.map2genome}"

rule reads_map2trans:
	input:
		'clean/{sample}.clean.fa'
	output:
		unmap2genome = 'mapping/unmap2genome/{sample}.unmap2genome.fa',
		map2trans = 'mapping/map2trans/{sample}.map2trans'
	params:
		DNA = config['DNA'],
		RNA = config['RNA']
	shell:
		 "bowtie -p5 -f -v1 {params.DNA} {input} --un {output.unmap2genome} > /dev/null; bowtie -p5 -f -v1 {params.RNA} {output.unmap2genome} |sed -r 's/\|.+\|\t/\t/' > {output.map2trans}"	

rule reads_num:
	input:
		fa = 'clean/{sample}.clean.fa',
		map2genome = 'mapping/map2genome/{sample}.map2genome',
		map2trans = 'mapping/map2trans/{sample}.map2trans'
	output:
		fa = 'stat/{sample}.clean.fa.num',
		map2genome = 'stat/{sample}.map2genome.num',
		map2trans = 'stat/{sample}.map2trans.num',
	shell:
		"echo -e \"{input.fa}\\t`grep -c '>' {input.fa}`\" > {output.fa}; echo -e \"{input.map2genome}\\t`cut -f1 {input.map2genome}|sort|uniq|wc -l`\" > {output.map2genome}; echo -e \"{input.map2trans}\\t`cut -f1 {input.map2trans}|sort|uniq|wc -l`\" > {output.map2trans}"

rule reads_stat:
	input:
		['stat/{sample}.clean.fa.num'.format(sample=x) for x in config['samples']],
		['stat/{sample}.map2genome.num'.format(sample=x) for x in config['samples']],
		['stat/{sample}.map2trans.num'.format(sample=x) for x in config['samples']]
	output:
		'stat/reads_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		"cat stat/*.num|sort > stat/all_num; {params.Rscript} bin/reads_stat.R"

rule get_iclip_site_genome:
	input:
		'mapping/map2genome/{sample}.map2genome'
	output:
		'sites/bed/{sample}.genome.bed'
	shell:
		"cat {input}|bin/upstream_1bp_pos_genome.py|sortBed > {output}"

rule get_iclip_site_trans:
	input:
		'mapping/map2trans/{sample}.map2trans'
	output:
		trans = 'sites/bed/{sample}.trans.bed',
		trans2genome = 'sites/bed/{sample}.trans2genome.bed'
	params:
		exon2genome = config['exon2genome']
	shell:
		"cat {input}|bin/upstream_1bp_pos_trans.py|sortBed > {output.trans}; bedtools intersect -a {output.trans} -b {params.exon2genome} -wa -wb|bin/trans_coor_convert.py|sortBed|uniq > {output.trans2genome}"

rule get_iclip_site_all:
	input:
		genome = 'sites/bed/{sample}.genome.bed',
		trans2genome = 'sites/bed/{sample}.trans2genome.bed'
	output:
		bed = 'sites/{sample}.sites.bed',
		bedGraph = 'track/{sample}.sites.bedGraph',
	shell:
		"cat {input.genome} {input.trans2genome}|sort -k4|groupBy -g 4 -c 5 -o count -full|awk -F '\\t' -vOFS='\\t' '{{print $1,$2,$3,$4,$8,$6}}'|sortBed > {output.bed}; cut -f1,2,3,5 {output.bed} > {output.bedGraph}"

rule bedgraph2tdf:
	input:
		'track/{sample}.sites.bedGraph'
	output:
		'track/{sample}.tdf'
	params:
		IGV = config['IGV']
	shell:
		"igvtools toTDF {input} {output} {params.IGV}"

rule reproduced_sites:
	input:
		['sites/{sample}.sites.bed'.format(sample=x) for x in config['samples']]
	output:
		tsv = 'reproduce/reproduce.counts.tsv',
		pdf = 'reproduce/reproduce.counts.pdf'
	shell:
		'bin/reproduced_sites.sh'

rule tag_clustering:
	input:
		['sites/{sample}.sites.bed'.format(sample=x) for x in config['samples']]
	output:
		'clusters/stat/tag_site_clusters.counts.tsv',
	params:
		feature = config['feature'],
		gene = config['GTF']
	shell:
		'bin/tag_clustering.sh {params.feature} {params.gene}'

