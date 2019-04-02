#!/bin/bash  
#$ -l h_rt=2:00:00
#$ -cwd 
#$ -j y
#$ -l h_vmem=6g
#$ -N compute_ld_scores_geneWindow
#$ -t 1-22

source /broad/software/scripts/useuse
reuse -q Python-2.7
reuse -q .anaconda-5.0.1

datetime=$(date)
echo Job $JOB_ID started $datetime >> 2_compute_ld_scores_geneWindow.log

window_size_in_kb=
label=
chrom=$SGE_TASK_ID

echo $chrom

python /cvar/jhlab/rebecca/ldsc/ldsc.py \
	--l2 \
	--bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
	--ld-wind-cm 1 \
	--annot ${label}_${window_size_in_kb}kb_geneWindow.${chrom}.annot.gz \
	--out ${label}_${window_size_in_kb}kb_geneWindow.${chrom} \
	--print-snps print_snps.txt

