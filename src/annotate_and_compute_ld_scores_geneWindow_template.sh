#!/bin/bash  
#$ -l h_rt=2:00:00
#$ -cwd 
#$ -j y
#$ -l h_vmem=6g
#$ -N annotate_and_compute_ld_scores_geneWindow
#$ -t 1-22

source /broad/software/scripts/useuse
reuse -q Python-2.7
reuse -q .anaconda2-5.0.1

chrom=$SGE_TASK_ID
label=depict_and_gtex_gene_intersection
window_size_in_kb=50
gene_boundary_path=../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt
bfile=ldsc/data/phase_3/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}
baseline_path=ldsc/data/phase_3/baseline_v1.1/baseline.
print_snps_path=ldsc/data/phase_3/1000G_EUR_Phase3_baseline/print_snps.txt

datetime=$(date)
echo Job $JOB_ID started $datetime 

python ../src/annotate_geneWindow.py \
        --label ${label} \
        --window_size_in_kb ${window_size_in_kb} \
        --curr_chromosome ${chrom} \
        --gene_boundary_file ${gene_boundary_path} \
        --baseline_file ${baseline_path}

python ldsc/ldsc.py \
        --l2 \
        --bfile ${bfile} \
        --ld-wind-cm 1 \
        --annot ${label}_${window_size_in_kb}kb_geneWindow.${chrom}.annot.gz \
        --out ${label}_${window_size_in_kb}kb_geneWindow.${chrom} \
        --print-snps ${print_snps_path}
