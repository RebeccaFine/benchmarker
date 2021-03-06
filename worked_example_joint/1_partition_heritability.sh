#!/bin/bash  
#$ -l h_rt=2:00:00
#$ -cwd 
#$ -j y
#$ -l h_vmem=8g
#$ -N partition_heritability
#$ -t 1
#$ -l os=RedHat7

source /broad/software/scripts/useuse
reuse -q Python-2.7
reuse -q .anaconda2-5.3.1

datetime=$(date)
echo Job $JOB_ID started $datetime >> partition_heritability.log

joblist=gwas_trait_list.txt

gwas_path=$(awk -F "\t" "NR==$SGE_TASK_ID+1 {print \$1}" ${joblist} )
trait=$(awk -F "\t" "NR==$SGE_TASK_ID+1 {print \$2}" ${joblist} )
label1=depict
label2=magma
input_directory_1=../worked_example_depict/
input_directory_2=../worked_example_magma/
output_label=depict_v_magma_joint
window_size_in_kb=50
output_directory=./
baseline_path=../../ldsc/data/phase_3/baseline_v1.1/baseline.
gene_window_path=../gene_window_files/depict_and_gtex_gene_intersection_50kb_geneWindow.
weights_path=../../ldsc/data/phase_3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frq_path=../../ldsc/data/phase_3/1000G_Phase3_frq/1000G.EUR.QC.
ldsc_script_path=../../ldsc/ldsc.py
normalize_script_path=../src/add_stats_to_results_files_joint.py # notice this one is different for the joint model

echo Job $JOB_ID started $datetime
echo job list: $joblist
echo gwas: $gwas_path
echo trait: $trait
echo label1: $label1
echo label2: $label2
echo input directory 1: $input_directory_1
echo input directory 2: $input_directory_2
echo output_label: $output_label
echo window size: $window_size_in_kb
echo gene window path: $gene_window_path
echo output directory: $output_directory
echo baseline: $baseline_path
echo weights: $weights_path
echo frq: $frq_path
echo ldsc_script_path: ${ldsc_script_path} 
echo normalize_script_path: ${normalize_script_path}

python ${ldsc_script_path} \
        --h2 ${gwas_path} \
        --ref-ld-chr ${baseline_path},${gene_window_path},${input_directory_1}/${trait}_${label1}_${window_size_in_kb}kb.,${input_directory_2}/${trait}_${label2}_${window_size_in_kb}kb. \
        --w-ld-chr ${weights_path} \
        --overlap-annot \
        --frqfile-chr ${frq_path} \
        --print-coefficients \
        --print-delete-vals \
        --out ${output_directory}/${trait}_${output_label}_${window_size_in_kb}kb

python ${normalize_script_path} \
	--trait ${trait} \
	--label1 ${label1} \
	--label2 ${label2} \
	--full_file_stem ${trait}_${output_label}_${window_size_in_kb}kb \
	--baseline_file_path ${baseline_path} \
	--input_file_directory ${output_directory} \
	--output_file_directory ${output_directory} 
