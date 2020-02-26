python 0_get_genes_from_ranked_list.py \
	--input_file /cvar/jhlab/rebecca/depict_runs/standard_5e-8/body_BMIz/body_BMIz_standard_5e-8_noChr@_geneprioritization_outside_input_loci_corrected.txt.gz \
	--gene_col_name 'Ensembl gene ID' \
	--sort_direction ascending \
	--col_to_sort_on 'P value' \
	--gene_boundary_file ../data/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt \
	--percentage_cutoff 10 \
	--trait body_BMIz \
	--label standard
