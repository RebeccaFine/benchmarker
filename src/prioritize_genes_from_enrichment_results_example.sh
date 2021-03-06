python ../src/0_prioritize_genes_from_enrichment_results.py \
	--output_directory . \
	--trait ibd_liu \
	--results_file /cvar/jhlab/rebecca/magma_runs/full_gene_sets/ibd_liu_noChr@_full_gene_sets.gcov.out \
	--gene_boundary_file ../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt \
	--set_definitions_file ../data/Top50GenesPerGeneSet_DEPICTGenes_Ensembl_depictAndGtexGeneIntersection.txt \
	--set_definitions_label Top50GenesPerGeneSet \
	--output_label test \
	--percentage_cutoff 10 \
	--results_file_set_col COVAR \
	--results_file_col_to_sort_on P \
	--sort_direction ascending \
	--results_file_separator '\s+'
