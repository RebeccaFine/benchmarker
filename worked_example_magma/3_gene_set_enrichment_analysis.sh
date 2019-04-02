#!/bin/bash  
#$ -cwd 
#$ -j y
#$ -l h_vmem=6g
#$ -N magma_genesets
#$ -l h_rt=01:00:00
#$ -l os=RedHat6
#$ -t 1-22

datetime=$(date)
echo Job $JOB_ID started $datetime >> magma_genesets.log

source /broad/software/scripts/useuse
reuse -q .anaconda-5.0.1

trait=bmi_magma
chrom=${SGE_TASK_ID}

echo $trait, $chrom, $gene_set_condition

../../magma/magma \
        --gene-results ${trait}_noChr${chrom}.genes.raw \
        --gene-covar ../data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z_GTExGenesOnly.txt onesided=greater \
        --out ${trait}_noChr${chrom}
