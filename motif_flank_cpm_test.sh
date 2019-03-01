source activate gmatic
GENOME=$HOME/Gmatic7/iCount/homo_sapiens.88.fa
#INPUT=MeCP2_KO_SE_RPKM_2
INPUT=MeCP2_KO_RBFOX2_KI_SE_RPKM_2_olp

#chr10:+:101232253:101232434:101229610:101229814:101237938:101238237

grep -v 'exonStart_0base' ${INPUT}|tr ':' '\t'|awk '{print $1"\t"$6"\t"$3"\t.\t.\t"$2"\n"$1"\t"$4"\t"$7"\t.\t.\t"$2}'|sed 's/chr//' > ${INPUT}_intron.bed
grep -v 'exonStart_0base' ${INPUT}|tr ':' '\t'|awk '{print $1"\t"$3-100"\t"$3"\t.\t.\t"$2"\n"$1"\t"$4"\t"$4+100"\t.\t.\t"$2}'|sed 's/chr//' > ${INPUT}_100nt.bed
grep -v 'exonStart_0base' ${INPUT}|tr ':' '\t'|awk '{print $1"\t"$3-300"\t"$3"\t.\t.\t"$2"\n"$1"\t"$4"\t"$4+300"\t.\t.\t"$2}'|sed 's/chr//' > ${INPUT}_300nt.bed

bedtools getfasta -fi $GENOME -bed ${INPUT}_intron.bed -s -tab -fo -|egrep 'GCATG|TGCAT' > ${INPUT}_intron_with_motif.tab
bedtools getfasta -fi $GENOME -bed ${INPUT}_100nt.bed -s -tab -fo -|egrep 'GCATG|TGCAT' > ${INPUT}_100nt_with_motif.tab
bedtools getfasta -fi $GENOME -bed ${INPUT}_300nt.bed -s -tab -fo -|egrep 'GCATG|TGCAT' > ${INPUT}_300nt_with_motif.tab

python find_motif.py ${INPUT}_intron_with_motif.tab 10 > ${INPUT}_intron_motif_flank10.bed
python find_motif.py ${INPUT}_intron_with_motif.tab 5 > ${INPUT}_intron_motif_flank5.bed

python find_motif.py ${INPUT}_100nt_with_motif.tab 10 > ${INPUT}_100nt_motif_flank10.bed
python find_motif.py ${INPUT}_100nt_with_motif.tab 5 > ${INPUT}_100nt_motif_flank5.bed

python find_motif.py ${INPUT}_300nt_with_motif.tab 10 > ${INPUT}_300nt_motif_flank10.bed
python find_motif.py ${INPUT}_300nt_with_motif.tab 5 > ${INPUT}_300nt_motif_flank5.bed

bedtools intersect -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank10_KO_rep1
bedtools intersect -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank10_KO_rep2
bedtools intersect -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank10_KO_rep3
bedtools intersect -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank10_WT_rep1
bedtools intersect -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank10_WT_rep2
bedtools intersect -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank10_WT_rep3

bedtools intersect -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank5_KO_rep1
bedtools intersect -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank5_KO_rep2
bedtools intersect -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank5_KO_rep3
bedtools intersect -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank5_WT_rep1
bedtools intersect -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank5_WT_rep2
bedtools intersect -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph -b ${INPUT}_intron_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_intron_motif_flank5_WT_rep3

bedtools intersect -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank10_KO_rep1
bedtools intersect -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank10_KO_rep2
bedtools intersect -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank10_KO_rep3
bedtools intersect -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank10_WT_rep1
bedtools intersect -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank10_WT_rep2
bedtools intersect -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank10_WT_rep3

bedtools intersect -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank5_KO_rep1
bedtools intersect -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank5_KO_rep2
bedtools intersect -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank5_KO_rep3
bedtools intersect -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank5_WT_rep1
bedtools intersect -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank5_WT_rep2
bedtools intersect -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph -b ${INPUT}_100nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_100nt_motif_flank5_WT_rep3

bedtools intersect -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank10_KO_rep1
bedtools intersect -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank10_KO_rep2
bedtools intersect -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank10_KO_rep3
bedtools intersect -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank10_WT_rep1
bedtools intersect -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank10_WT_rep2
bedtools intersect -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank10.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank10_WT_rep3

bedtools intersect -a peaks/MeCP2_KO_rep1.TGGTCA_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank5_KO_rep1
bedtools intersect -a peaks/MeCP2_KO_rep2.CACTGT_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank5_KO_rep2
bedtools intersect -a peaks/MeCP2_KO_rep3.ATTGGC_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank5_KO_rep3
bedtools intersect -a peaks/MeCP2_WT_rep1.CGTGAT_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank5_WT_rep1
bedtools intersect -a peaks/MeCP2_WT_rep2.ACATCG_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank5_WT_rep2
bedtools intersect -a peaks/MeCP2_WT_rep3.GCCTAA_reads_unique_peaks.bedGraph -b ${INPUT}_300nt_motif_flank5.bed -wa -wb|awk '{print $8"\t"$4}' > ${INPUT}_300nt_motif_flank5_WT_rep3

Rscript motif_flank_cpm_t-test.R $INPUT
rm ${INPUT}_*nt*
