configfile: "config.yaml"

rule all:
	input:
		expand('fastq/{sample}.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_fastqc.html', sample=config['samples']),
		expand('trim/{sample}.trim3adapter.fastq.gz', sample=config['samples']),
		expand('fastqc/trim/{sample}.trim3adapter_fastqc.html', sample=config['samples']),
		expand('dedup/{sample}.dedup.fa', sample=config['samples']),
		expand('mapping/{sample}', sample=config['samples']),
		expand('mapping/{sample}/Aligned.sortedByCoord.out.bam', sample=config['samples']),
		expand('xlsites/{sample}_reads_unique.bed', sample=config['samples']),
		expand('xlsites/{sample}_reads_multiple.bed', sample=config['samples']),
		expand('xlsites/{sample}_reads_skipped.bam', sample=config['samples']),
		expand('peaks/{sample}_reads_unique_peaks.bed', sample=config['samples']),
		expand('peaks/{sample}_reads_unique_peaks_scores.tsv', sample=config['samples']),
		'reproduce/crosslink_sites_reproduce.pdf',
		'merge/merge_ALL_crosslink_sites.bed',
		'merge/merge_ALL_crosslink_sites_sig.bed',
		'merge/merge_ALL_crosslink_sites_sig_clusters.bed',
		'merge/merge_ALL_crosslink_sites_anno_biotype.tab',
		'merge/merge_ALL_crosslink_sites_anno_geneid.tab',
		'merge/merge_ALL_crosslink_sites_sig_anno_biotype.tab',
		'merge/merge_ALL_crosslink_sites_sig_anno_geneid.tab',
		'merge/merge_ALL_crosslink_sites_summary.tab',
		'merge/merge_ALL_crosslink_sites_sig_summary.tab',
		'figure/crosslink_sites_distr.pdf',
		'table/stats_tab.tsv',
		'track/merge_ALL_crosslink_sites.tdf',
		'track/merge_ALL_crosslink_sites_sig.tdf',


rule split_fastq:
	input:
		config['raw_reads']
	output:
		'fastq/{sample}.fastq.gz'
	shell:
		'bin/split_fastq_by_barcode.py {input} {output}'

rule fastqc_raw:
	input:
		'fastq/{sample}.fastq.gz'
	output:
		'fastqc/raw/{sample}_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -o fastqc/raw {input}'

rule trim3adapter:
	input:
		'fastq/{sample}.fastq.gz'
	output:
		'trim/{sample}.trim3adapter.fastq.gz'
	params:
		conda = config['conda_path'],
		adapt = config['adapt']
	shell:
		 "{params.conda}/cutadapt -f fastq -a {params.adapt} -m15 {input} | gzip -c > {output}"

rule fastqc_trim:
	input:
		'trim/{sample}.trim3adapter.fastq.gz'
	output:
		'fastqc/trim/{sample}.trim3adapter_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -o fastqc/trim {input}'

rule remove_PCR_duplicates:
	input:
		'trim/{sample}.trim3adapter.fastq.gz'
	output:
		'dedup/{sample}.dedup.fa'
	params:
		conda = config['conda_path']
	shell:
		 "zcat {input} | {params.conda}/fastq_to_fasta -Q33 -r -n | bin/split_barcode.py | sort | uniq | cut -f1 | bin/seq2fasta.py > {output}"

rule mapping:
	input:
		'dedup/{sample}.dedup.fa'
	output:
		dir = 'mapping/{sample}',
		bam = 'mapping/{sample}/Aligned.sortedByCoord.out.bam'
	params:
		conda = config['conda_path'],
		genome = config['genome'],
		annotation = config['annotation'],
		mismatch = config['mismatch'],
		cpu = config['cpu']
	shell:
		'mkdir -p {output.dir} && export ICOUNT_TMP_ROOT={output.dir}.tmp && {params.conda}/iCount mapstar -mis {params.mismatch} --threads {params.cpu} {input} {params.genome} {output.dir} --annotation {params.annotation}'

rule xlsites_reads:
	input:
		'mapping/{sample}/Aligned.sortedByCoord.out.bam'
	output:
		'xlsites/{sample}_reads_unique.bed',
		'xlsites/{sample}_reads_multiple.bed',
		'xlsites/{sample}_reads_skipped.bam'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/iCount xlsites {input} {output} --group_by start --quant reads'

rule significant_xlsites:
	input:
		'xlsites/{sample}_reads_unique.bed',
	output:
		bed = 'peaks/{sample}_reads_unique_peaks.bed',
		tsv = 'peaks/{sample}_reads_unique_peaks_scores.tsv'
	params:
		conda = config['conda_path'],
		segmentation = config['segmentation'],
		FDR = config['FDR']
	shell:
		'{params.conda}/iCount peaks {params.segmentation} {input} {output.bed} --scores {output.tsv} --fdr {params.FDR}'

rule reproduce:
	input:
		["peaks/{sample}_reads_unique_peaks.bed".format(sample=x) for x in config['samples']]
	output:
		'reproduce/crosslink_sites_reproduce.pdf'
	shell:
		'bin/reproduced_sites.sh'

rule merge_xlsites:
	input:
		["xlsites/{sample}_reads_unique.bed".format(sample=x) for x in config['samples']]
	output:
		'merge/merge_ALL_crosslink_sites.bed'
	shell:
		'bin/merge_rep_sites.sh'

rule significant_merged_xlsites:
	input:
		'merge/merge_ALL_crosslink_sites.bed'
	output:
		'merge/merge_ALL_crosslink_sites_sig.bed'
	params:
		conda = config['conda_path'],
		segmentation = config['segmentation'],
		FDR = config['FDR']
	shell:
		'{params.conda}/iCount peaks {params.segmentation} {input} {output.bed} --fdr {params.FDR}'

rule clusters:
	input:
		'merge/merge_ALL_crosslink_sites_sig.bed'
	output:
		'merge/merge_ALL_crosslink_sites_sig_clusters.bed'
	params:
		conda = config['conda_path'],
		cluster_dist = config['cluster_dist']
	shell:
		'{params.conda}/bedtools cluster -d {params.cluster_dist} -i {input} > {output}'
 
rule xlsites_sig_annotate_biotype:
	input:
		'merge/merge_ALL_crosslink_sites_sig.bed'
	output:
		'merge/merge_ALL_crosslink_sites_sig_anno_biotype.tab',
	params:
		conda = config['conda_path'],
		segmentation = config['segmentation']
	shell:
		'{params.conda}/iCount annotate --subtype biotype {params.segmentation} {input} {output}'

rule xlsites_sig_annotate_gene_id:
	input:
		'merge/merge_ALL_crosslink_sites_sig.bed'
	output:
		'merge/merge_ALL_crosslink_sites_sig_anno_geneid.tab',
	params:
		conda = config['conda_path'],
		segmentation = config['segmentation']
	shell:
		'{params.conda}/iCount annotate --subtype gene_id {params.segmentation} {input} {output}'

rule xlsites_ALL_annotate_biotype:
	input:
		'merge/merge_ALL_crosslink_sites.bed'
	output:
		'merge/merge_ALL_crosslink_sites_anno_biotype.tab',
	params:
		conda = config['conda_path'],
		segmentation = config['segmentation']
	shell:
		'{params.conda}/iCount annotate --subtype biotype {params.segmentation} {input} {output}'

rule xlsites_ALL_annotate_gene_id:
	input:
		'merge/merge_ALL_crosslink_sites.bed'
	output:
		'merge/merge_ALL_crosslink_sites_anno_geneid.tab',
	params:
		conda = config['conda_path'],
		segmentation = config['segmentation']
	shell:
		'{params.conda}/iCount annotate --subtype gene_id {params.segmentation} {input} {output}'

rule xlsites_summary_all:
	input:
		'merge/merge_ALL_crosslink_sites_anno_geneid.tab',
	output:
		'merge/merge_ALL_crosslink_sites_summary.tab',
	params:
		conda = config['conda_path'],
	shell:
		'bin/sites_summary.sh {input} {output}'

rule xlsites_summary_sig:
	input:
		'merge/merge_ALL_crosslink_sites_sig_anno_geneid.tab',
	output:
		sig = 'merge/merge_ALL_crosslink_sites_sig_summary.tab',
	params:
		conda = config['conda_path'],
	shell:
		'bin/sites_summary.sh {input} {output}'

rule xlsites_distr:
	input:
		'merge/merge_ALL_crosslink_sites_sig_summary.tab',
		'merge/merge_ALL_crosslink_sites_summary.tab'
	output:
		'figure/crosslink_sites_distr.pdf'
	shell:
		'bin/xlsites_distr.R {input} {output}'

rule stats_tab:
	input:
		["xlsites/{sample}_reads_unique.bed".format(sample=x) for x in config['samples']],
	output:
		'table/stats_tab.tsv'		
	shell:
		'bin/stats.sh'

rule bed2bedgraph:
	input:
		sig = 'merge/merge_ALL_crosslink_sites_sig.bed',
		ALL = 'merge/merge_ALL_crosslink_sites.bed'
	output:
		sig = 'track/merge_ALL_crosslink_sites_sig.bedGraph',
		ALL = 'track/merge_ALL_crosslink_sites.bedGraph'
	shell:
		"cut -f1,2,3,5 {input.sig} > {output.sig}; cut -f1,2,3,5 {input.ALL} > {output.ALL}"

rule bedgraph2tdf:
	input:
		sig = 'track/merge_ALL_crosslink_sites_sig.bedGraph',
		ALL = 'track/merge_ALL_crosslink_sites.bedGraph'
	output:
		sig = 'track/merge_ALL_crosslink_sites_sig.tdf',
		ALL = 'track/merge_ALL_crosslink_sites.tdf'
	params:
		IGV = config['IGV']
	shell:
		"igvtools toTDF {input.sig} {output.sig} {params.IGV}; igvtools toTDF {input.ALL} {output.ALL} {params.IGV}"

