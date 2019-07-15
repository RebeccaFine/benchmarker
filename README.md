# Benchmarker
Author: Rebecca S. Fine

Note 7/15/19: The code for generating annotations + LD scores for gene windows has been simplified. The generate_gene_window/ folder has been removed, and the necessary script and template (annotate_geneWindow.py and annotate_and_compute_ld_scores_geneWindow_template.sh) are in the src/ directory.

This repository contains a collection of scripts people may find useful if they wish to implement our Benchmarker strategy, described in this manuscript: https://www.cell.com/ajhg/fulltext/S0002-9297(19)30146-6. If you use these, please cite:  
Fine et al. (2019) American Journal of Human Genetics

If you have questions, you can email me at rebeccafine@g.harvard.edu or open an issue in the repository.

Software:   
* LDSC software: https://github.com/bulik/ldsc (version 1.0) - we also recommend looking over the wiki here
* DEPICT software: https://data.broadinstitute.org/mpg/depict/ (release 194) (see important note about DEPICT below)
* MAGMA software: https://ctg.cncr.nl/software/magma (version 1.06b)  


Dependencies:  
Pandas (>= 0.17)  
NumPy  
SciPy  
intervaltree  
argparse  

Data:  
* Gene boundary annotation file used in the manuscript (based on genes found in both DEPICT and GTEx): data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt   
* LD scores and other files for LDSC: https://data.broadinstitute.org/alkesgroup/LDSCORE (more information at https://github.com/bulik/ldsc)
* 1000G_Phase3_baseline_v1.1_ldscores.tgz 
* 1000G_Phase3_frq.tgz  
* 1000G_Phase3_plinkfiles.tgz  
* 1000G_Phase3_weights_hm3_no_MHC.tgz 
* FTP site: ftp://ftp.broadinstitute.org/outgoing/benchmarker_data/
    * Binarized gene sets (*_DEPICTGenes_Ensembl_depictAndGtexGeneIntersection.txt)
    * Gene window files (gene_window_depict_and_gtex_gene_intersection_50kb.tar.gz)
    * MAGMA .annot file (depictAndGtexGeneIntersection_noMHC.gene.loc, depictAndGtexGeneIntersection_noMHC.gene.annot)

We have provided several use-case examples:

1) General implementation for gene prioritization results 
2) DEPICT
3) MAGMA (this is also a fairly general implementation for enrichment-based results)

For all examples, download https://data.broadinstitute.org/alkesgroup/UKBB/body_BMIz.sumstats.gz

# Important Note about DEPICT
There are some versions of DEPICT in which the _geneprioritization_outside_input_loci.txt file (which is the relevant one for prioritization) is missing a tab between the 'P value' and 'False discovery rate' columns for genes with FDR >= 0.20. Please check your output files and fix this if this problem is present.

# Important first step for any of these analyses: calculating LD scores for gene windows
## EDITED for simplicity 7/15/19

The first step is to generate what we refer to as "gene window" LD scores. Specifically, we will annotate all SNPs within 50 kb of any gene in your defined gene boundary as "1" and everything else as "0". This will serve as a control when we partition heritability (i.e. controlling for the fact that SNPs near/within any gene might be enriched for heritability). A template for this is in src/annotate_and_compute_ld_scores_geneWindow_template.sh. The necessary paths and parameters are listed at the top and can be changed for your particular situation. (The script calls src/annotate_geneWindow.py to perform the actual annotation.)
	 
The necessary files for running these can be downloaded from https://data.broadinstitute.org/alkesgroup/LDSCORE/. For more information, check out the LDSC wiki: https://github.com/bulik/ldsc/. In this case, we use gene boundaries in GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt (provided on the Github and in the FTP site); these are the gene boundaries used in the AJHG paper, which represent genes in both DEPICT and GTEx excluding the MHC.

Note that you can decide to use a window size other than 50 kb by changing the parameter in the 1_annotate_geneWindow.sh script. If you do that, just make sure that that window size is consistent for *both* your gene window annotation and your prioritization annotation.

Also note that the gene window files used in the manuscript are on the FTP site (gene_window_depict_and_gtex_gene_intersection_50kb.tar.gz). These are based on the input files called in the template script.

# General use case (gene prioritization)

For this, we assume you have already run your method of interest in a leave-one-chromosome-out fashion.  For each of these runs, prioritize the top 10% of genes on the withheld chromosome and put them into a single list of gene IDs.

Once you have generated your list of prioritized input genes, you will need to 1) annotate all SNPs in the baseline model as to whether they are prioritized or not, 2) calculate LD scores for the annotation you've created, 3) partition the heritability for your GWAS of interest, and 4) normalize the partitioned heritability by average per-SNP heritability for the trait. 

I have structured this code so it is easy to test multiple GWASes at once, which I hope will be useful if you are testing performance of a particular method.  I have also included a script, generate_submission_scripts.py, which takes as input a single config file (generate_submission_scripts_example) and will produce 3 shell scripts that cover all of the other steps in this section.

Note that template submission files are in a format appropriate for submitting a job as a task array on an UGER-based system. This should be relatively straightforward to convert to any other type of queueing system (e.g. LSF).

You will need to process your input GWAS for LDSC with their munge_sumstats.py script (see their wiki for details); this will format your GWAS data with the appropriate input columns. Here is an example command for munging our BMI summary statistics.

```
python ldsc/munge_sumstats.py \
        --sumstats body_BMIz.sumstats.gz \
        --out body_BMIz_forLDSC.sumstats.gz \
        --merge-alleles ldsc/data/w_hm3.snplist
```

Another technical note to be aware of is that for the GWAS analyzed with LMM methods (e.g. BOLT-LMM), the effective N is larger than the actual N. This will not matter for any values that are normalized by overall trait h2 (e.g. normalized tau, proportion of h2 explained), but it does affect raw tau values (and of course their standard error). Note that in our manuscript, we focus exclusively on normalized tau, and we recommend others do the same; however, in the interest of thoroughness, we mention this here.

If you wish to analyze these values, you will need to divide them by (Neff/N).  For the GWAS in our manuscript from BOLT-LMM, these values can be found in the README of https://data.broadinstitute.org/alkesgroup/UKBB. For a more detailed explanation of this issue, see the supplement of "Functional architecture of low-frequency variants highlights strength of negative selection across coding and non-coding annotations" (Gazal et al. 2018).

**A full example is in benchmarker/worked_example_from_gene_list/. Go to this folder to execute the following steps.**

## Necessary input files
* gwas_trait_list.txt: A (tab-delimited) file with the following headers: gwas_path, trait_name. Here, list the path to each GWAS you will use for LDSC (processed appropriately with munge_stats.py from the LDSC pipeline, see above) and a label for the trait.
```
gwas_path       trait_name
../example_files/body_BMIz_forLDSC.sumstats.gz bmi
../example_files/ibd_liu_forLDSC.sumstats.gz    ibd
```
* list of priortiized genes (in this case, examples are in bmi_genes.txt and ibd_genes.txt)

## Automated pipeline
    python src/generate_submission_scripts.py generate_submission_scripts_example.cfg

This script takes as input a single cfg file (src/generate_submission_scripts_example.cfg). The cfg file contains all the arguments necessary to create 3 shell scripts, listed in the next section.

Most of the arguments in this cfg file should be self-explanatory, but here are a few notes:
* The convention in this pipeline is that all output files take the form [TRAIT]\_[LABEL]\_[WINDOW_SIZE_IN_KB]kb. This naming structure makes it easy to parallelize testing many GWASes at once.
* Your input list of genes may or may not have a header. Define this as a True or False in 'prio_gene_list_has_header = '. If there is a header, specify the column name containing the gene names in 'gene_col_name'
* Make sure you have defined a control annotation for the gene window (see "Calculating LD scores for gene windows) above -- this is what goes in 'gene_window_path = ../gene_window_files/depict_and_gtex_gene_intersection_50kb_geneWindow.'
* All LDSC files can be downloaded from the https://data.broadinstitute.org/alkesgroup/LDSCORE/ ftp site (and more information about LDSC on their github page at https://github.com/bulik/ldsc)
* You almost certainly will not need to change the paths under the '# these will generally remain unchanged' section - these are just the names of the templates the script is using to create the task arrays and our in-house scripts for ${annotate_script_path} and ${normalize_script_path} (src/annotate_input_gene_list.py and src/add_stats_to_results_files.py, respectively)


### Output scripts from generate_submission_scripts.py:
**Note that the expected output scripts can be found in worked_example_from_gene_list/expected_output_scripts/ (for purposes of comparison/inspection if you're trying to generate these yourself)**

1) [LABEL]\_annotateAndComputeLDScores.sh: annotating SNPs + calculating LD scores. 
* Format is UGER task array (run as qsub [LABEL]\_annotateAndComputeLDScores.sh)
* Outputs: .annot.gz files, .l2.ldscore.gz files (these are standard LDSC outputs)
* ${ldsc_script_path} is the path to ldsc.py

```
# File contents:
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
```

2) [LABEL]\_partitioningHeritability: partitioning heritability + tau normalization. 
* Format is UGER task array (run as qsub [LABEL]\_partitioningHeritability)
* Outputs: .results, .results_withStats files. 
* The .results file is the original LDSC output file for partitioning h2
* The .results_withStats file has additional columns for normalized tau/tau standard error
* Note that LDSC also produces .part_delete files; these are used to calculate the jackknifed standard errors for the partitioned h2 and will be important if you are interested in generating a p-value for comparing methods directly

```
# File contents:
python ${ldsc_script_path} \
--h2 ${gwas_path} \
--ref-ld-chr ${baseline_path},${gene_window_path},${output_directory}/${trait}_${label}_${window_size_in_kb}kb. \
--w-ld-chr ${weights_path} \
--overlap-annot \
--frqfile-chr ${frq_path} \
--print-coefficients \
--print-delete-vals \
--out ${output_directory}/${trait}_${label}_${window_size_in_kb}kb

python ${normalize_script_path} \
--trait ${trait} \
--label ${label} \
--full_file_stem ${trait}_${label}_${window_size_in_kb}kb \
--baseline_file_path ${baseline_path} \
--input_file_directory ${output_directory} \
--output_file_directory ${output_directory}
```

3) [LABEL]\_concatenate.sh  concatenating .results files into a single file. Format is a script that can be run directly on the command line. (run as ./[LABEL]\_concatenate.sh)
* Output: .concatenatedResults file. This combines the results from each GWAS into a single output file (convenient if you'd like to compare multiple GWASes)

### Summary of steps:
```
python ../src/generate_submission_scripts.py generate_submission_scripts_example.cfg
qsub test_annotateAndComputeLDScores.sh
qsub test_partitionHeritability.sh
./test_concatenate.sh 
```

# DEPICT use case
**Example in benchmarker/worked_example_depict/.**

We assume you have already run DEPICT on your GWAS of interest in a leave-one-chromosome-out fashion (this example is based on body_BMIz.sumstats.gz).  Example cfgs for this can be found in worked_example_depict/body_BMIz_standard_5e-8_noChr*.cfg (you will need to create 22 cfgs, in which each one specifies withholding a different chromosome). The critical piece is the following argument in the cfg:
```
# Chromosome to be left-out, leave empty if all chromosomes should be included (for bencharmking purposes)
leave_out_chr: 1
```

Some DEPICT prioritization output files are provided here for the purposes of this example: body_BMIz_standard_5e-8_noChr*_geneprioritization_outside_input_loci_corrected.txt.gz

Steps:


1) ./1_get_genes_from_ranked_list_example.sh

python ../src/0_get_genes_from_ranked_list.py \
    --input_file ../example_files/body_BMIz_standard_5e-8_noChr@_geneprioritization_outside_input_loci_corrected.txt.gz \
    --gene_col_name 'Ensembl gene ID' \
    --sort_direction ascending \
    --col_to_sort_on 'P value' \
    --gene_boundary_file ../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt \
    --percentage_cutoff 10 \
    --trait bmi \
    --label depict_gene_sets

 This code goes through each provided prioritization file, sorts by p-value, and takes the top 10% of genes for each chromosome.  You can modify the percentage of genes desired by changing '--percentage cutoff'.  The '--trait' and '--label' flags will define the resulting outputs, which take the form [trait]_[label]_prioTop[prioritized_percentage]Genes.txt.  This output file contains a list genes to be prioritized (no header).

2) After this, the steps are the same as the general use case; a relevant cfg file is provided in this folder (2_generate_submission_scripts_depict.cfg)

**Note that the expected output scripts can be found in worked_example_depict/expected_output_scripts/ (for purposes of comparison/inspection if you're trying to generate these yourself)**

```
	python ../src/generate_submission_scripts.py 2_generate_submission_scripts_depict.cfg 
```
Generates:
```
depict_annotateAndComputeLDScores.sh # generates LD score files
depict_partitionHeritability.sh # partitions heritability and normalizes tau 
depict_concatenate.sh # groups multiple tested GWAS into a single file (if there are any)
```

### Summary of steps:
```
./1_get_genes_from_ranked_list_example.sh
python ../src/generate_submission_scripts.py 2_generate_submission_scripts_depict.cfg 
qsub depict_annotateAndComputeLDScores.sh
qsub depict_partitionHeritability.sh
./depict_concatenate.sh 
```
# MAGMA use case

Prior to this, make sure GWAS of interest is in MAGMA-friendly format (columns = SNP, P, N, where SNP is rsID).

You will also need a .annot file, which lists all genes and the SNPs within their boundaries.
To calculate the particular .annot files we used (depictAndGtexGeneIntersection_noMHC.genes.annot):
* Download data/depictAndGtexGeneIntersection_noMHC.gene.loc (from this repository)
* Download MAGMA and assciated data (1kg/g1000*)
 * The 1kg/g1000* files can be downloaded from the "Reference data" table; we use European data, but obviously you should select the population that matches your data
* Run ../magma --annotate --snp-loc ../1kg/g1000_eur.bim --gene-loc depictAndGtexGeneIntersection_noMHC.gene.loc --out depictAndGtexGeneIntersection_noMHC

You can also use their provided gene.loc file, which contains more genes (NCBI37.3.genes.loc; this is found in their "Auxiliary files" table, and you will likely want the build 37 data to match the build of LDSC), or make your own based on your preferred list of gene boundaries.  From their or any other .gene.loc file, you can create a .annot file the same way as described above (with the --annotate --snp-loc command referenced earlier). Detailed instructions for this are given in the MAGMA manual.


1) Calculate MAGMA gene-level results from GWAS
	qsub 1_gene_level_results.sh 
	
	```
	# File contents:
	../../magma/magma \
    --bfile ../../magma/1kg/g1000_eur \
    --pval ../example_files/body_BMIz_forMAGMA.txt  ncol=N --gene-annot ../data/depictAndGtexGeneIntersection_noMHC.genes.annot --out bmi_magma
    ```
* Produces bmi_magma.genes.raw, bmi_magma.genes.out (gene-level summary statistic information)

2) Filter each gene-level result to remove a chromosome (argument: label of .genes.raw output files)
```
	python 2_withhold_chromosome.py bmi_magma
```
* Note that for this, we edit the .genes.raw data (rather than the .genes.out data), as this is what is used in the gene set enrichment step

3) Perform MAGMA gene set enrichment analysis
	qsub 3_gene_set_enrichment_analysis.sh
	
	```
	# File contents:
	../../magma/magma \
        --gene-results ${trait}_noChr${chrom}.genes.raw \
        --gene-covar ../data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z_GTExGenesOnly.txt onesided=greater \
        --out ${trait}_noChr${chrom}
    ```

4) From enrichment results, prioritize the top genes on the withheld chromosome
	./4_prioritize_genes_from_enrichment_results_example.sh

	```
	# File contents:
	python ../src/prioritize_genes_from_enrichment_results.py \
    --output_directory . \
    --trait ${trait} \
    --results_file ${trait}_magma_noChr@.gcov.out \
    --gene_boundary_file ../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt \
    --set_definitions_file ../data/Top50GenesPerGeneSet_DEPICTGenes_Ensembl_depictAndGtexGeneIntersection.txt \
    --set_definitions_label Top50GenesPerGeneSet \
    --output_label magma \
    --percentage_cutoff 10 \
    --results_file_set_col COVAR \
    --results_file_col_to_sort_on P \ 
    --sort_direction ascending \
    --results_file_separator '\s+'
    ```
* This produces a file of prioritized genes (with a header), including not only the gene IDs but the chromosome and the enriched gene set the prioritized gene came from. Note that this script can also be used for converting other types of enrichment results to gene prioritization -- see below for details.

5) After this, the steps are the same as the general use case; a relevant cfg file is provided in this folder (5_generate_submission_scripts_example_magma.cfg):
	python ../src/generate_submission_scripts.py 5_generate_submission_scripts_example.cfg


**Note that the expected output scripts can be found in worked_example_magma/expected_output_scripts/ (for purposes of comparison/inspection if you're trying to generate these yourself)**
	Then submit:
```
qsub magma_annotateAndComputeLDScores.sh 
qsub magma_partitionHeritability.sh
./magma_concatenate.sh
```

Summary of steps:
```
qsub 1_gene_level_results.sh 
python 2_withhold_chromosome.py bmi_magma
qsub 3_gene_set_enrichment_analysis.sh
./4_prioritize_genes_from_enrichment_results_example.sh
python ../src/generate_submission_scripts.py 5_generate_submission_scripts_example.cfg
qsub magma_annotateAndComputeLDScores.sh 
qsub magma_partitionHeritability.sh
./magma_concatenate.sh

```
# Making a joint model and calculating p-values between methods

If you'd like to compare two methods directly, we recommend putting them together in a joint model. 

**Example in benchmarker/worked_example_joint/**

## Partitioning h2 in a joint model

Modeling two methods jointly is quite straightforward.  After calculating the LD scores for each individually (from earlier steps), you can then partition the heritability by using both sets of LD scores in the '--ref-ld-chr' flag (commma-separated).  In this example, label1 and label2 are the labels differentiating the methods (e.g. "depict" and "magma").  All arguments besides --ref-ld-chr and --normalize_script_path are the same as for the separate models (don't forget you still need the gene window annotation!).

```
qsub 1_partition_heritability.sh 

# File contents:
python ${ldsc_script_path} \
        --h2 ${gwas_path} \
        --ref-ld-chr ${baseline_path},${gene_window_path},${input_directory_1}/${trait}_${label1}_${window_size_in_kb}kb.,${input_directory_2}/${trait}_${label2}_${window_size_in_kb}kb. \
        --w-ld-chr ${weights_path} \
        --overlap-annot \
        --frqfile-chr ${frq_path} \
        --print-coefficients \
        --print-delete-vals \
        --out ${output_directory}/${trait}_${output_label}_${window_size_in_kb}kb
```

We also use a different script to normalize the tau values (add_stats_to_results_files_joint.py rather than add_stats_to_results_files.py). This can be called as follows:

```
qsub 2_depict_v_magma_concatenate.sh

# File contents:
normalize_script_path=../src/add_stats_to_results_files_joint.py
python ../src/${normalize_script_path} \
        --trait ${trait} \
        --label1 ${label1} \
        --label2 ${label2} \
        --full_file_stem ${trait}_${output_label}_${window_size_in_kb}kb \
        --baseline_file_path ${baseline_path} \
        --input_file_directory ${output_directory} \
        --output_file_directory ${output_directory}
``` 

## Calculating p-values
Calculating p-values within a joint model is a simple extension of the jackknife procedure LDSC uses to calculate its standard errors.  Instead of calculating a standard error around a mean, we're going to calculate the standard error around the difference between two means (i.e. tau_1 - tau_2).  To do this, we need to use the .part_delete outputs from LDSC, which contain the block jackknife estimates from the original standard error calculation. We provide src/compare_pvals_joint.py, which performs this calculation based on the .concatenatedResults file.  (For this reason, there is a column in the .concatenatedResults file containing the path to the .part_delete file.) 

IMPORTANT NOTE: this script is written assuming that there are exactly four arguments to ref-ld-chr in the following order: baseline_model, gene_window, ldScoresFromMethod1, ldScoresFromMethod2.  If this assumption is violated, results will not be accurate -- the script has to know which columns of the .part_delete file to use, and it is specifically assuming that the annotations of interest are columns 54 and 55 of this file (which will be true if the assumption is met).  For each trait, the output contains a line with "Theta / theta_hatJ"; this should generally be close to 1.  If it isn't, check that your files are in the expected order.

```
./3_compare_pvals_joint.sh

# File contents:
python ../src/compare_pvals_joint.py \
        --results_directory . \
        --file_prefix depict_v_magma_joint \
        --output_dir . \
        --output_label depict_v_magma_joint
```

Summary of steps:
```
qsub 1_partition_heritability.sh 
./2_depict_v_magma_concatenate.sh
3_compare_pvals_joint.sh
```


# A note about enrichment --> prioritization

```
../src/prioritize_genes_from_enrichment_results.py
```

Used to prioritize genes based on withhold-one-chromosome runs of your enrichment method of interest (we already saw a usage case for MAGMA).

This script works for enrichment results rather than gene prioritization results -- it will convert enrichment results to gene prioritization. This is the process that was used to analyze MAGMA in our manuscript. This requires some extra inputs relative to the previous script, specifically defining gene members of the sets that were used in the analysis.

The input must be 22 files where the left-out chromosome is in the name of the file. In the results_file flag, replace the chromosome number with @.

The gene set definition file (--set_definitions_file) is formatted as follows. Each line represents one gene set (defined as whatever enrichment you are performing -- e.g. could be sets of genes expressed in different tissues(. First column is the name of the set.  Member genes are listed after the name of the set (tab-delimited).

SET_1   ENS1    ENS2    ENS3...  
SET_2   ENS3    ENS5    ENS3...  

Note that at the moment, there is some implicit ranking in the gene set definitions, where genes listed earlier are considered better than those listed later. This only matters in the case where the script is approaching the 10% cutoff for number of genes -- if it can only prioritize only some genes in the gene set before hitting the cutoff, it will always go from first-listed to last.

```
        python 0_prioritize_genes_from_enrichment_results.py \
        --results_file
        --set_definitions_file
        --percentage_cutoff # we recommend 10, but you can change this
        --results_file_set_col # name of the column with gene set IDs in the results file
        --results_file_col_to_sort_on # p-value, FDR, rank, etc.
        --sort_direction # "ascending" or "descending"
        --results_file_separator #e.g. '\t' '\s+' ' ' 
        --gene_boundary_file ../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt # or gene boundary file of your choice
        --output_label
        --trait
        --output_directory
```



