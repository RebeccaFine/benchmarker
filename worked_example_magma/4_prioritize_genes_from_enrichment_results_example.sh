trait=bmi
# note that full output label will be TRAIT_OUTPUTLABEL

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
