import pandas as pd
import numpy as np
import scipy.stats
import argparse
import datetime
parser = argparse.ArgumentParser()
parser.add_argument('--file_prefix', help = 'Prefix for .results file for the first concatenated se of results to look at')
parser.add_argument('--results_directory', help = 'Directory containing the first set of results')
parser.add_argument('--output_dir')
parser.add_argument('--output_label')

args = parser.parse_args()

now = datetime.datetime.now()
print 'Started:', now.strftime("%Y-%m-%d %H:%M")
print '\n'

print 'Arguments:'
print 'reading in IDs and positions...'
for arg in vars(args):
    print arg, ':', getattr(args, arg)
print '\n'

print 'Assuming the prioritization annotation column for both sets of results is L2_2, which is the last column of the block jackknife file.'
file_prefix = args.file_prefix
results_directory = args.results_directory
if not results_directory.endswith('/'):
    results_directory = results_directory + '/' 
output_dir = args.output_dir
output_label = args.output_label

def block_jackknife( merged_df ):
    jackknife_file = results_directory + '/' + merged_df[ 'jackknife_file_path_' + label_1 ]

    jackknife= pd.read_csv( jackknife_file, sep = ' ', header = None)
    
    #  normalize jackknife coefficient values by the same normalizaion factor as used previously (i.e. h2/M, average per-SNP h2)
    normalization_factor_1 = merged_df[ 'normalization_factor_' + label_1 ]
    normalization_factor_2 = merged_df[ 'normalization_factor_' + label_2 ]
    
    j1_normalized = jackknife.iloc[:,54] * normalization_factor_1 
    j2_normalized = jackknife.iloc[:,55] * normalization_factor_2
    
    diffs = j1_normalized - j2_normalized
    print '\n'

    obs_coefficient_1 = float( merged_df[ 'tau_normalized_' + label_1] )
    obs_coefficient_2 = float(  merged_df[ 'tau_normalized_' + label_2] )
    obs_difference = obs_coefficient_1 - obs_coefficient_2

    theta = obs_difference
    theta_j = np.array(diffs)
    theta_hatJ= np.sum(theta-theta_j)+sum(theta_j/200.0)
    tau = 200.0 * theta-199* theta_j
    print 'Trait:', merged_df['Trait']
    print 'theta (= observed difference):', theta
    print 'theta_hatJ:', theta_hatJ
    standard_error = np.sqrt( (1/200.0)* sum((tau-theta_hatJ)**2/(199.0)))
    print 'Theta / theta_hatJ:', theta / theta_hatJ

    if obs_difference > 0:
        merged_df['Higher_Value'] = label_1
    else:
        merged_df['Higher_Value'] = label_2

    z = obs_difference / standard_error
    print 'standard error:', standard_error
    print 'Z-score:', z
    print 'P-value (two-tailed):', scipy.stats.norm.sf(abs(z))*2 #twosided
    #scipy.stats.norm.sf(abs(z))*2

    merged_df['Observed_difference'] = obs_difference
    merged_df['Jackknife_std_error'] = standard_error
    merged_df['Jackknife_zscore'] = z
    merged_df['Jackknife_pvalue'] = scipy.stats.norm.sf(abs(z))*2 

    return merged_df

# read in results files 1 and 2 - these should be concatenated files where each row is one trait's prioritization result

results_file = '{results_directory}/{file_prefix}.concatenatedResults'.format( results_directory = results_directory, file_prefix = file_prefix)

print 'results file:', results_file
results = pd.read_csv( results_file, sep = '\t')

x = set(results[ results.Category == 'L2_2']['Label'].tolist())
y = set(results[ results.Category == 'L2_3']['Label'].tolist())
if len(x) > 1:
    raise Exception('Labels inconsistent for L2_2')
if len(y) > 1:
    raise Exception('Labels inconsistent for L2_3')

label_1 = list(x)[0]
label_2 = list(y)[0]
print 'inferred label 1:', label_1
print 'inferred label 2:', label_2
results_reshaped = results.pivot(index='Trait', columns='Label')
results_reshaped.columns = results_reshaped.columns.map('_'.join)  
results_reshaped.reset_index(inplace = True)
print results_reshaped.head(5)


results_reshaped = results_reshaped.apply( block_jackknife, axis = 1 )

outname = output_dir + '/' + output_label + '_PVals.txt'
results_reshaped.to_csv(outname, sep = '\t', index = False)
print '\n'
print 'Output in', outname

