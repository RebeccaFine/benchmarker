import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--label')
parser.add_argument('--joblist')
parser.add_argument('--output_directory')

args = parser.parse_args()

label = args.label
joblist = args.joblist
output_directory = args.output_directory

outname = output_directory + '/' + label + '.concatenatedResults'

overall_df = pd.DataFrame()
with open(joblist) as x:
    x.readline()
    for line in x:
        gwas, trait = line.strip().split('\t')
        filename =  '{output_directory}/{trait}_{label}_50kb.results_withStats'.format( trait = trait, label = label, output_directory = output_directory)
        df = pd.read_csv(filename , sep = '\t')
        overall_df = pd.concat( [overall_df, df])

overall_df.to_csv( outname, index = False, sep = '\t')
