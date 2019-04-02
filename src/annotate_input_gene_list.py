import pandas as pd
import gzip
import sys
import numpy as np
#sys.path.append( '/broad/software/free/Linux/redhat_6_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/lib/python2.7/site-packages' )  
import intervaltree as it
import argparse
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--prio_gene_list', help = 'List of prioritized genes. Can be either a single list of gene IDs with no header or a file with multiple columns with headers, in which one column contains gene IDs. In the latter case, specify the name of the column containing gene IDs with --gene_col_name.')
parser.add_argument('--curr_chromosome',  help = 'Chromosome to prioritize on')
parser.add_argument('--gene_boundary_file', help = 'File with all gene IDs, chromosomes, and start/stop positions. Gene ID header should be ensembl_id, chromosome header should be ensembl_chromosome, ensembl_bp_start, ensembl_bp_end (tab-delimited). May contain other columns')
parser.add_argument('--window_size_in_kb', help = 'Number of kilobases on either side of each gene to include for SNP assignment')
parser.add_argument('--output_label', help = 'Label for output. Output file format is $trait_$label_$window_size_in_kb')
parser.add_argument('--baseline_file', help = 'Path to baseline model for LDSC')
parser.add_argument('--prio_gene_list_has_header', help = 'True or False; indicates whether there is a header row in the gene prioritization file')
parser.add_argument('--gene_col_name', help ='If there is a header row, name of the column containing gene ID') 

args = parser.parse_args()

prio_gene_list = args.prio_gene_list
curr_chromosome = args.curr_chromosome
gene_boundary_file = args.gene_boundary_file
window_size_in_kb = int(args.window_size_in_kb)
output_label = args.output_label
baseline_file = args.baseline_file
prio_gene_list_has_header = args.prio_gene_list_has_header
if prio_gene_list_has_header.lower() == 'false':
    prio_gene_list_has_header_bool = False
elif prio_gene_list_has_header.lower() == 'true':
    prio_gene_list_has_header_bool = True
else:
    raise Exception('prio_gene_list_has_header_bool should be True or False')
gene_col_name = args.gene_col_name

now = datetime.datetime.now()
print 'Started:', now.strftime("%Y-%m-%d %H:%M")
print '\n'
"""
print 'Input gene list:', prio_gene_list
print 'Chromosome:', curr_chromosome
print 'Gene boundary file:', gene_boundary_file # must have columns ensembl_id, ensembl_bp_start, ensembl_bp_end, ensembl_chromosome
print 'Window size (in kb):', window_size_in_kb
print 'Output label:', output_label
print 'prio_gene_list_has_header?', prio_gene_list_has_header_bool
if  prio_gene_list_has_header_bool:
    print 'prio_gene_list column header:', gene_col_name
print '\n'
"""

# https://stackoverflow.com/questions/34992524/print-command-line-arguments-with-argparse
print 'Arguments:'
print 'reading in IDs and positions...'
for arg in vars(args):
    print arg, ':', getattr(args, arg)

#IDsAndPositions = {} # key = Ensembl ID, value = (chrom, start, end)
gene_boundaries = pd.read_csv( gene_boundary_file, sep = '\t', dtype = {'ensembl_chromosome' : np.int32, 'ensembl_bp_start': np.int32,  'ensembl_bp_end' : np.int32 })
# restrict to current chromosome
gene_boundaries = gene_boundaries[ gene_boundaries.ensembl_chromosome == int(curr_chromosome) ]

gene_boundaries['tup'] = zip(gene_boundaries['ensembl_chromosome'], gene_boundaries['ensembl_bp_start'], gene_boundaries['ensembl_bp_end'])
IDsAndPositions = gene_boundaries.set_index('ensembl_id')['tup'].to_dict()


# make list of prioritized genes on relevant chromosome

if prio_gene_list_has_header_bool == False:
    input_prio_genes = pd.read_csv( prio_gene_list, sep= '\t', header = None).iloc[:,0].tolist()  # okay if all chromosomes are included - Ensembl dict will filter out anything not on current chromosome
else:
    input_prio_genes = pd.read_csv( prio_gene_list, sep= '\t')[ gene_col_name ].tolist()

print '\n'
print 'length of input gene list (may include genes on all chromosomes):', len(input_prio_genes)
prio_genes = [ gene for gene in input_prio_genes if gene in IDsAndPositions ] 
print 'number of genes on chromosome %s: %s' %(str(curr_chromosome), str(len(prio_genes)))
print 'percentage of genes prioritzed on chromosome {chrom} = ({len_prio} / {num_genes}) = {fraction}'.format( len_prio = len(prio_genes), num_genes = len(IDsAndPositions), fraction = len(prio_genes)/float(len(IDsAndPositions)), chrom = curr_chromosome)

gene_output = '{output_label}_chrom{chrom}.txt'.format( output_label = output_label, chrom = curr_chromosome )
print 'Prioritized genes for given chromosome in', gene_output
with open(gene_output, 'w') as x:
    for gene in prio_genes:
        x.write(gene + '\n' )

# make interval trees of prioritized genes
prio_tree = it.IntervalTree()
for gene in prio_genes:
    chrom, start, stop = IDsAndPositions[ gene ]
    start_withWindow = int(start) - (window_size_in_kb*1000)
    stop_withWindow = int(stop) + (window_size_in_kb*1000)
    prio_tree[ start_withWindow: stop_withWindow] = gene

print '\n'

# go through each SNP in annotation file and annotate as prioritized or not
not_found = set()
output_annot_file = '{output_label}_{window_size_in_kb}kb.{chrom}.annot.gz'.format( output_label = output_label, window_size_in_kb = window_size_in_kb, chrom = curr_chromosome )

with gzip.open( output_annot_file , 'wb') as output_annot:
    with gzip.open( baseline_file ) as annot_file:
        annot_file.readline()
        
        cutoff_label='prio'
        header = ['CHR','BP','SNP', 'CM', cutoff_label ]

        output_annot.write( '\t'.join(header) + '\n')
        
        for line in annot_file:
            chrom, pos, rsid, cm, base = line.strip().split('\t')[0:5]
            
            if prio_tree.search( int(pos) ):
                prio_annot = '1'
            else:
                prio_annot = '0'
            
            output_annot.write( '\t'.join([chrom, pos, rsid, cm, prio_annot]) + '\n')            
print '\n'
print 'Output in:', output_annot_file
