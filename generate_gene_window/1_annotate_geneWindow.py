import pandas as pd
import gzip
import intervaltree as it
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--curr_chromosome')
parser.add_argument('--label')
parser.add_argument('--window_size_in_kb')
parser.add_argument('--gene_boundary_file') # /cvar/jhlab/rebecca/depict_versions/depict_download_10-3-16/DEPICT/data/mapping_and_annotation_files/GPL570ProbeENSGInfo+HGNC_reformatted_noMHC_depictAndGtexGeneIntersection_RF.txt

args = parser.parse_args()

curr_chromosome = args.curr_chromosome
label = args.label
window_size_in_kb = int( args.window_size_in_kb )
gene_boundary_file = args.gene_boundary_file

print 'Label:', label
print '\n'
print 'Window size:', window_size_in_kb
print'\n'
print 'Chromosome:', curr_chromosome
print '\n'
print 'Gene boundary file:', gene_boundary_file

IDsAndPositions = {} # key = Ensembl ID, value = (chrom, start, end, strand)

print 'reading in IDs and positions...'
gene_boundaries = pd.read_csv( gene_boundary_file, sep = '\t', dtype = {'ensembl_chromosome' : np.int32, 'ensembl_bp_start': np.int32,  'ensembl_bp_end': np.int32 })
gene_boundaries = gene_boundaries[ gene_boundaries.ensembl_chromosome == int(curr_chromosome) ]

gene_boundaries['tup'] = zip(gene_boundaries['ensembl_chromosome'], gene_boundaries['ensembl_bp_start'], gene_boundaries['ensembl_bp_end'])
IDsAndPositions = gene_boundaries.set_index('ensembl_id')['tup'].to_dict()
print 'number of genes:', len(IDsAndPositions)

# make interval tree of all genes on relevant chromosome
gene_tree = it.IntervalTree()
for gene in IDsAndPositions:
    start = IDsAndPositions[gene][1] - (window_size_in_kb*1000) # include window size on either side of interval
    end = IDsAndPositions[gene][2] + (window_size_in_kb*1000)
    gene_tree[ start:end ] = gene # add 1 because upper range of interval is exclusive
#print 'all:\n'
#print gene_tree
print '\n'
print 'tree length, all genes:', len(gene_tree)
print '\n'

# go through each SNP in annotation file and annotate as prioritized or nonprioritized 
not_found = set()
output_annot_file = '%s_%skb_geneWindow.%s.annot.gz' %(label, str(window_size_in_kb), str(curr_chromosome))
original_baseline_file = '/cvar/jhlab/rebecca/ldsc/data/phase_3/baseline_v1.1/baseline.%s.annot.gz' %str(curr_chromosome)

with gzip.open( output_annot_file, 'wb') as output_annot:
    with gzip.open( original_baseline_file ) as annot_file:
        annot_file.readline()
        
        header = ['CHR','BP','SNP', 'CM', 'in_any_gene_window']

        output_annot.write( '\t'.join(header) + '\n')
        
        for line in annot_file:
            chrom, pos, rsid, cm, base = line.strip().split('\t')[0:5]
            if gene_tree.search( int(pos) ):
                gene_annot = '1'
            else:
                gene_annot = '0'
            
            output_annot.write( '\t'.join([chrom, pos, rsid, cm, gene_annot]) + '\n')            
print '\n'
print 'Output in:', output_annot_file
