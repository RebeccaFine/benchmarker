#!/bin/bash
#$ -l h_rt=2:00:00
#$ -cwd
#$ -j y
#$ -l h_vmem=6g
#$ -N annotate_and_calc_ld_scores
#$ -t 1-44
#$ -l os=RedHat7

source /broad/software/scripts/useuse
reuse -q Python-2.7
reuse -q .anaconda2-5.3.1

datetime=$(date)
echo Job $JOB_ID started $datetime >> annotate_and_calc_ld_scores.log

joblist=test_jobList.txt

# from listed job parameters
trait=$(awk -F "\t" "NR==$SGE_TASK_ID+1 {print \$2}" ${joblist} )
curr_chromosome=$(awk -F "\t" "NR==$SGE_TASK_ID+1 {print \$3}" ${joblist} )

# consistent
prio_gene_list=@_genes.txt
output_directory=./
label=test
window_size_in_kb=50
baseline_path=../../ldsc/data/phase_3/baseline_v1.1/baseline.${curr_chromosome}.
gene_boundary_path=../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt
bfile_path=../../ldsc/data/phase_3/1000G_EUR_Phase3_plink/1000G.EUR.QC.${curr_chromosome}
print_snps_path=../../ldsc/data/phase_3/1000G_EUR_Phase3_baseline/print_snps.txt
prio_gene_list_has_header=False
gene_col_name=NA
ldsc_script_path=/cvar/jhlab/rebecca/ldsc/ldsc.py
annotate_script_path=../src/annotate_input_gene_list.py

echo Job $JOB_ID started $datetime
echo trait:, $trait
echo label:, $label
input_gene_list=$(echo ${prio_gene_list} | sed "s#@#${trait}#g")
echo input gene list: $input_gene_list
echo curr_chromosome: $curr_chromosome
echo output directory:, $output_directory
echo window size: $window_size_in_kb
echo prio_gene_list_has_header: $prio_gene_list_has_header
echo gene_col_name: $gene_col_name

echo bfile path: $bfile_path
echo baseline path: $baseline_path
echo print snps path: $print_snps_path
echo ldsc_script_path: $ldsc_script_path
echo annotate_script_path: $annotate_script_path
echo '\n\n'

python ${annotate_script_path} \
--prio_gene_list ${input_gene_list} \
--curr_chromosome ${curr_chromosome} \
--gene_boundary_file  ${gene_boundary_path} \
--output_label ${output_directory}/${trait}_${label} \
--window_size_in_kb 50 \
--baseline_file ${baseline_path}annot.gz \
--prio_gene_list_has_header ${prio_gene_list_has_header} \
--gene_col_name ${gene_col_name}

python ${ldsc_script_path} \
--l2 \
--bfile ${bfile_path} \
--ld-wind-cm 1 \
--annot ${output_directory}/${trait}_${label}_${window_size_in_kb}kb.${curr_chromosome}.annot.gz \
--out ${output_directory}/${trait}_${label}_${window_size_in_kb}kb.${curr_chromosome} \
--print-snps ${print_snps_path}
