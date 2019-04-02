import pandas as pd
import sys, gzip
import argparse
import numpy as np
import datetime

now = datetime.datetime.now()
print 'Started:', now.strftime("%Y-%m-%d %H:%M")
print '\n'
print '\n'

parser = argparse.ArgumentParser()
parser.add_argument('--trait', help = 'Trait of interest. This will be the first part of the output label for all output files')
parser.add_argument('--results_file', help = 'Enrichment result to convert to prioritized genes')
parser.add_argument('--gene_boundary_file', help = 'File with all gene IDs in both Ensembl and HGNC, plus chromosomes. Gene ID header should be ensembl_id, chromosome header should be ensembl_chromosome, HGNC header should be HGNC (tab-delimited). May contain other columns. Note that number of genes to prioritize per chromosome will be defined based on this file. Therefore, make sure these boundaries are defined only for the genes that could potentially be prioritized from your data set.')
parser.add_argument('--set_definitions_file', help = 'File where each line represents one gene set (defined as whatever enrichment you are performing -- e.g. could be sets of genes expressed in different tissues. First column is the name of the set.  Member genes are listed after the name of the set (tab-delimited)') #'/cvar/jhlab/rebecca/magma/gene_sets/depict_and_gtex_intersection/%s_DEPICTGenes_Ensembl_depictAndGtexGeneIntersection.txt
parser.add_argument('--set_definitions_label', help = 'Label for the set definitions used')
parser.add_argument('--output_label', help = 'Note that final output format is {trait}_{output_label}_{set_definitions_label}_Top{percentage_cutoff}PercentGenes_Genes')
parser.add_argument('--percentage_cutoff', help = "Percentage of genes on each chromosome to prioritiz. Note this is defined based on the number of genes per chromosome from the gene boundary file")
parser.add_argument('--results_file_set_col', help = 'Column in the results file containing the names of the enriched gene sets')
parser.add_argument('--results_file_col_to_sort_on',help = "Column in enrichment file to sort gene sets on (e.g. p-values)")
parser.add_argument('--sort_direction', help = 'Direction to sort gene sets in prioritization file (ascending or descending). P-values would be ascending, for example')
parser.add_argument('--results_file_separator', help = 'Column separator for results file. Can be \t, ' ', or \s+')
parser.add_argument('--output_directory')

args = parser.parse_args()
print 'Arguments:'
print 'reading in IDs and positions...'
for arg in vars(args):
    print arg, ':', getattr(args, arg)
print '\n'


output_directory = args.output_directory
if not output_directory.endswith('/'):
    output_directory = output_directory + '/'
trait = args.trait
results_file = args.results_file
gene_boundary_file = args.gene_boundary_file
output_label = args.output_label
set_definitions_file = args.set_definitions_file
set_definitions_label = args.set_definitions_label
percentage_cutoff  = int(args.percentage_cutoff )
results_file_col_to_sort_on = args.results_file_col_to_sort_on
sort_direction = args.sort_direction
results_file_separator = args.results_file_separator 
results_file_set_col = args.results_file_set_col

if results_file_separator not in ['\t',' ', '\s+']:
    raise Exception('--results_file_separator must be \\t, "' '", or \s+')

if sort_direction == 'ascending':
    sort_direction_bool = True
elif sort_direction == 'descending':
    sort_direction_bool = False
else:
    raise Exception('sort_direction must be "ascending" or "descending"')

gene_boundaries = pd.read_csv( gene_boundary_file, sep = '\t', dtype = {'ensembl_chromosome' : np.int32, 'ensembl_bp_start': np.int32,  'ensembl_bp_end' : np.int32 })

gene_boundaries = gene_boundaries[ gene_boundaries.ensembl_chromosome.isin( range(1,23) )].copy()
ensembl_to_hugo_dict = gene_boundaries.set_index('ensembl_id')[ 'HGNC' ].to_dict()
hugo_to_ensembl_dict = gene_boundaries.set_index('HGNC')[ 'ensembl_id' ].to_dict()
ensembl_dict = gene_boundaries.groupby('ensembl_chromosome')['ensembl_id'].apply(lambda x: x.tolist()).to_dict()


# define genes in each set
gene_set_dict = {}
with open( set_definitions_file ) as x:
    for line in x:
        gene_set, genes = line.strip().split('\t')[0], line.strip().split('\t')[1:]
        gene_set_dict[ gene_set ] = genes


output_df = pd.DataFrame()
print trait
for chrom in range(1,23):
    print 'Chromosome:', chrom
    filename = results_file.replace('@',str(chrom) )

    results = pd.read_csv( filename, sep = results_file_separator, comment = '#')
    #print results.head(5)
    results.sort_values( results_file_col_to_sort_on, inplace = True, ascending = sort_direction_bool )
    results.reset_index(inplace =True, drop = True )
    ordered_genesets = results[ results_file_set_col ].tolist()
    num_genes = int(round( len(ensembl_dict[ int(chrom) ]) * float(percentage_cutoff)/100.0))
    print 'Number of genes on this chromosome defined in gene boundary file:', len(ensembl_dict[ int(chrom) ])
    print 'Number of genes to prioritize:', num_genes
    
    final_gene_list = []
    final_gene_set_list = []
        
    for gene_set in ordered_genesets:
        if len(final_gene_list) > num_genes:
                break

        if gene_set in gene_set_dict: # if it's not in the dict, that means it didn't have enough genes passing the filter
            genes = gene_set_dict[ gene_set ]
            for gene in genes:
                if gene not in final_gene_list and gene in  ensembl_dict[ int(chrom) ]:
                    final_gene_list.append( gene )
                    final_gene_set_list.append( gene_set )
        else:
            print 'gene set not found (presumably failed filter):', gene_set

    final_gene_list = final_gene_list[0: num_genes ]
    final_gene_set_list = final_gene_set_list[0: num_genes]
    chromosome_list = [chrom] * num_genes
    trait_list = [trait] * num_genes
    label_list = [ set_definitions_label ] * num_genes

    output_df = output_df.append( pd.DataFrame( [trait_list, label_list, chromosome_list, final_gene_list, final_gene_set_list]).T)

output_df.columns = ['Trait','Label','Chromosome','Ensembl_Gene','Set_Name']

output_df['HUGO_Gene_tmp'] = output_df.Ensembl_Gene.map( ensembl_to_hugo_dict )

def fix_hugo_gene(df):
    if df['HUGO_Gene_tmp'] == '-':
        return df['Ensembl_Gene']
    else:
        return df['HUGO_Gene_tmp']

output_df['HUGO_Gene'] = output_df.apply( fix_hugo_gene, 1)

column_order = ['Trait','Label','Chromosome','Ensembl_Gene','HUGO_Gene','Set_Name']
output_df = output_df[ column_order ]


output_df_label = '{output_directory}/{trait}_{output_label}_{set_definitions_label}_Top{percentage_cutoff}PercentGenes_Genes.txt'.format( chrom = str(chrom), trait = trait, output_label = output_label, percentage_cutoff = percentage_cutoff, set_definitions_label = set_definitions_label, output_directory = output_directory)

print 'Output in:', output_df_label 
output_df.to_csv(output_df_label, sep = '\t', index = False)

