import argparse
import pandas as pd
import datetime

# Reads in gene prioritization files from leave-one-chromosome-out runs and creates the top-10% list of genes for input to Benchmarker (taking the top 10% of genes on chromosome 1 from the leave-out-chromosome-1 run, the top 10% of genes on chromosome 2 from the leave-out-chromosome-2 run, etc.) 10% can be changed to desired percentage.


parser = argparse.ArgumentParser()
parser.add_argument( '--input_file', help = 'Prioritization files. This generally will be 22 separate files, one per chromosome (where the left-out chromosome is in the name of the file; replace the chromosome number with @ in the file name -- genes on the withheld chromosome will be used for prioritization). \
 Alternatively, you may specify a single input file that contains genes on all chromosomes *IF* the values in that file come from leave-one-out runs (i.e. the p-values for all genes listed for chromosome 1 were calculated in a run which left out chromosome 1, the p-values for all genes listed for chromosome 2 were calculated in a run which left out chromosome 2, etc.')
parser.add_argument('--gene_col_name', help = 'Name of column containing genes in prioritization file')
parser.add_argument('--col_to_sort_on', help = 'Column in prioritization file to sort genes on (e.g. p-values)')
parser.add_argument('--sort_direction', help = 'Direction to sort genes in prioritization file (ascending or descending). P-values would be ascending, for example')
parser.add_argument('--chrom_col_name', help = 'If chromosome column is present in prioritization file, put the name here. If absent, chromosome will be determined from the gene boundary file')
parser.add_argument('--gene_boundary_file', help = 'File with all gene IDs and chromosomes. Gene ID header should be ensembl_id, chromosome header should be ensembl_chromosome (tab-delimited). May contain other columns')
parser.add_argument('--percentage_cutoff', help = 'Percentage of genes to prioritize on each chromosome. We recommend 10')
parser.add_argument('--trait', help = 'Name of trait prioritization is done for. This will be part of the output file name')
parser.add_argument('--label', help = 'Label for output file. Output file format will be $trait_$label_Genes.txt')

args = parser.parse_args()
print 'Arguments:'
print 'reading in IDs and positions...'
for arg in vars(args):
    print arg, ':', getattr(args, arg)
print '\n'



input_file_name = args.input_file
sort_direction = args.sort_direction
col_to_sort_on = args.col_to_sort_on
gene_boundary_file = args.gene_boundary_file
percentage_cutoff = args.percentage_cutoff
trait = args.trait
label = args.label
gene_col_name = args.gene_col_name
chrom_col_name = args.chrom_col_name

now = datetime.datetime.now()
print 'Started:', now.strftime("%Y-%m-%d %H:%M")
print '\n'
print '\n'

if sort_direction == 'ascending':
    sort_direction_bool = True
elif sort_direction == 'descending':
    sort_direction_bool = False
else:
    raise Exception('sort_direction must be "ascending" or "descending"')

gene_boundaries = pd.read_csv( gene_boundary_file, sep = '\t')
all_genes = set(gene_boundaries['ensembl_id'].tolist())
print 'number of genes in boundary file:', len(all_genes)
genes_and_chrom_dict = gene_boundaries.set_index('ensembl_id')['ensembl_chromosome'].to_dict()

"""
# determine what number = relevant percentage
print '\n'
print 'Before removal:', prio_file.shape[0]
prio_file = prio_file[ prio_file[ 'Ensembl gene ID'].isin( IDsAndPositions ) ].copy()
print 'After removal:', prio_file.shape[0]
print 'Number of genes on current chromosome in given output file:', prio_file.shape[0]
"""

# if all chromosomes are in input file
if input_file_name.endswith('gz'):
    compression = 'gzip'
else:
    compression = None
overall_df = pd.DataFrame()
if '@' not in input_file_name:
    # read in input file
    input_file = pd.read_csv( input_file_name, sep = '\t', compression = compression)
    input_file = input_file[ input_file[ gene_col_name ].isin( all_genes ) ].copy()

    # if no chromosome column, add one
    if chrom_col_name == None:
        print 'Adding chromosome column'
        chrom_col_name = 'Chromosome'
        input_file[ 'Chromosome' ] = input_file[ gene_col_name ].map( genes_and_chrom_dict )
        input_file[ chrom_col_name ] = input_file[ chrom_col_name ].astype('int64')
        #print input_file.head(5)
   
    for chrom in range(1, 23):
        input_file_chrom = input_file[ input_file[ chrom_col_name] == int(chrom) ].copy()
        input_file_chrom.sort_values( col_to_sort_on, inplace = True, ascending = sort_direction_bool)
        input_file_chrom.reset_index(inplace = True, drop = True)
        num_rows = input_file_chrom.shape[0] 
        num_genes_to_use = int(round(float(num_rows) * (float(percentage_cutoff)/100.0 )))
        prio_genes = input_file_chrom.iloc[0: num_genes_to_use ][ gene_col_name  ]
        #prio_genes = input_file_chrom.iloc[0: num_genes_to_use ]['Ensembl gene ID'].tolist()
        print 'number of prioritized genes (top {percentage_cutoff}%): {num}'.format( percentage_cutoff =str(percentage_cutoff) , num = str(prio_genes.shape[0]) )
        overall_df = pd.concat( [overall_df, prio_genes ] )

# if each chromosome in a separate input file
else:
    for chrom in range(1,23):
        filename = input_file_name.replace('@', str(chrom) )
        input_file_chrom = pd.read_csv( filename, sep = '\t', compression = compression)
        input_file_chrom = input_file_chrom[ input_file_chrom[ gene_col_name ].isin( all_genes ) ].copy()
        
        # if no chromosome column, add one
        if chrom_col_name == None:
            print 'Adding chromosome column'
            chrom_col_name = 'Chromosome'
            input_file_chrom[ 'Chromosome' ] = input_file_chrom[ gene_col_name ].map( genes_and_chrom_dict )
            input_file_chrom[ chrom_col_name ] = input_file_chrom[ chrom_col_name ].astype('int64')
            #print input_file_chrom.head(5)

        input_file_chrom = input_file_chrom[ input_file_chrom[ chrom_col_name ]  == chrom ].copy()
        input_file_chrom.sort_values( col_to_sort_on, inplace = True, ascending = sort_direction_bool)
        input_file_chrom.reset_index(inplace = True, drop = True)
        num_rows = input_file_chrom.shape[0] 
        num_genes_to_use = int(round(float(num_rows) * (float(percentage_cutoff)/100.0 )))
        prio_genes = input_file_chrom.iloc[0: num_genes_to_use ][ gene_col_name ]
        #prio_genes = input_file_chrom.iloc[0: num_genes_to_use ]['Ensembl gene ID'].tolist()
        print 'number of prioritized genes (top {percentage_cutoff}%): {num}'.format( percentage_cutoff =str(percentage_cutoff) , num = str(prio_genes.shape[0]) )
        overall_df = pd.concat( [overall_df, prio_genes ] )

gene_output = '%s_%s_prioTop%sPercent_Genes.txt' %(trait, label, percentage_cutoff )
overall_df.to_csv( gene_output, sep = '\t', index = False, header = False)
print 'output in:', gene_output
print '\n'

