#!/bin/bash  
#$ -cwd 
#$ -j y
#$ -l h_vmem=10g
#$ -N magma_genes
#$ -l h_rt=1:00:00
#$ -l os=RedHat6

source /broad/software/scripts/useuse
use gcc-4.4.7
use Anaconda

datetime=$(date)
echo Job $JOB_ID started $datetime >> magma_genes.log


../../magma/magma \
    --bfile ../../magma/1kg/g1000_eur \
    --pval ../example_files/body_BMIz_forMAGMA.txt  ncol=N --gene-annot ../data/depictAndGtexGeneIntersection_noMHC.genes.annot --out bmi_magma
