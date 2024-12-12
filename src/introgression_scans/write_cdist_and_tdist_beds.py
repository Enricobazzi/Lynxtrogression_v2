import pandas as pd
import numpy as np
import os
import argparse

# parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Find distance from centromere and telomere')
    parser.add_argument('--ibed', type=str, help='Input bed file')
    parser.add_argument('--ct_file', type=str, help='Centromere file')
    parser.add_argument('--otfile', type=str, help='Output bed of telomere distances')
    parser.add_argument('--ocfile', type=str, help='Output bed of centromere distances')
    return parser.parse_args()

args = parse_args()
ibed = args.ibed
ct_file = args.ct_file
otfile = args.otfile
ocfile = args.ocfile

# read in the bed file
df = pd.read_csv(ibed, sep='\t', header=None, names=['chrom', 'start', 'end'],
                 dtype={'chrom': str, 'start': float, 'end': float})
# add the midpoint column
df['mid'] = (df['start'] + df['end']) / 2

# read in the centromere file
ct_data = pd.read_csv(ct_file, sep='\t', header=None, names=['chrom', 'length', 'centromere'],
                      dtype={'chrom': str, 'length': float, 'centromere': float})
# create dictionary of centromere and length
c_dict = dict(zip(ct_data['chrom'], ct_data['centromere']))
t_dict = dict(zip(ct_data['chrom'], ct_data['length']))

# add centromere and telomere columns to the bed file
df['centromere'] = df['chrom'].map(c_dict)
df['length'] = df['chrom'].map(t_dict)

# calculate centromere distance
df['c_dist'] = np.abs(df['mid'] - df['centromere'])
df['t_dist'] = np.min([df['mid'], df['length'] - df['mid']], axis=0)

# write out centromere and telomere distances bed files
df[['chrom', 'start', 'end', 'c_dist']].to_csv(ocfile, sep='\t', header=False, index=False)
df[['chrom', 'start', 'end', 't_dist']].to_csv(otfile, sep='\t', header=False, index=False)
