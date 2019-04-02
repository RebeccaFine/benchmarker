import pandas as pd
import numpy as np
import sys, os
import argparse
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--input_file_directory', help = 'Directory where .results files are located')
parser.add_argument('--output_file_directory', help = 'Directory where you would like the output results to go. Recommended to be the same as --input_file_directory, as it keeps this file in the same location as the jackknife files that will be used later for p-value calculation')
parser.add_argument('--full_file_stem', help = 'Full prefix for .results file of interest. Note that in the standard pipeline, this consists of TRAIT_LABEL_WINDOWSIZEkb; it is included as a flag for flexibility, but we recommend using this naming system.')
parser.add_argument('--trait', help = 'Name of trait being evaluated')
parser.add_argument('--label', help = 'Label for analysis (should not include trait name). In the suggested pipline, this is the same label as used throughout the LD score and h2 calculation')
parser.add_argument('--baseline_file_path', help = 'Path + prefix to LDSC baseline file used for partitioning h2 (this is only used to count the number of common SNPs included in the analysis)')

args = parser.parse_args()

trait = args.trait
label = args.label
output_file_directory = args.output_file_directory
if not output_file_directory.endswith('/'):
    output_file_directory = output_file_directory + '/'
input_file_directory = args.input_file_directory
if not input_file_directory.endswith('/'):
    input_file_directory = input_file_directory + '/'
full_file_stem = args.full_file_stem
baseline_file = args.baseline_file_path

now = datetime.datetime.now()
print 'Started:', now.strftime("%Y-%m-%d %H:%M")
print '\n'
print '\n'
# https://stackoverflow.com/questions/34992524/print-command-line-arguments-with-argparse
print 'Arguments:'
print 'reading in IDs and positions...'
for arg in vars(args):
    print arg, ':', getattr(args, arg)
print '\n'

# get list of SNPs with MAF > 5%
cmd = "cat %s*.l2.M_5_50  | cut -f1 | awk '{s+=$1} END {print s}'" %baseline_file
total_num_snps_included = float(os.popen(cmd).read().strip())

# get h2 estimate
results_file = input_file_directory + '/' + full_file_stem + '.results'
original_results =pd.read_csv( results_file, sep = '\t')

log_file = input_file_directory + '/' + full_file_stem + '.log'
with open(log_file) as x:
    for line in x:
        if line.startswith('Total Observed scale h2'):
            total_h2 = float(line.strip().split(':')[1].split(' ')[1])
            total_h2_std_error = float( line.strip().split('(')[1].strip(')'))

print 'total observed_scale h2:', total_h2

output_file = '{output_file_directory}/{full_file_stem}.results_withStats'.format( full_file_stem = full_file_stem , output_file_directory = output_file_directory)
original_results_file = pd.read_csv( results_file, sep = '\t')

# calculate normalized tau
original_results[ 'tau_normalized' ] = original_results['Coefficient'] * total_num_snps_included / float( total_h2 )
original_results['tau_standard_error_normalized']  = original_results['Coefficient_std_error'] * total_num_snps_included / float( total_h2 )
original_results['normalization_factor'] = float(total_num_snps_included) / float( total_h2 )
original_results['jackknife_file_path'] =  '{full_file_stem}.part_delete'.format( full_file_stem = full_file_stem)

# subset only to prioritized 
original_results = original_results[ original_results.Category == 'L2_2'].copy()
print 'Assuming L2_2 is the prioritized annotation'
original_columns = original_results.columns.tolist()
original_results['Trait'] = trait
original_results['Label'] = label
original_results = original_results[ ['Trait', 'Label'] + original_columns ]
original_results.to_csv( output_file, sep = '\t', index = False)
print 'Output in %s' %output_file
