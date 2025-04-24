import pandas as pd
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Write dmax values to bed file')
    parser.add_argument('--ivcf', help='Path to the phased vcf file with the haplotypes')
    parser.add_argument('--ibed', help='Path to the bed file with the windows')
    parser.add_argument('--pop_sizes', help='Comma-separated list of the number of individuals in each population')
    parser.add_argument('--obed', help='Path to the output bed file')
    return parser.parse_args()

args = parse_args()

vcf_file = args.ivcf
bed_file = args.ibed
pop_sizes = args.pop_sizes

def ham_dist(bin1, bin2):
    if bin1 == bin2:
        return 0
    else:
    # only works on binary strings (0 or 1) of equal length
        return bin(int(bin1, 2) ^ int(bin2, 2)).count('1')

nhaps = sum(map(int, pop_sizes.split(',')))

# read vcf file
vcf = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
# split the genotype columns into haplotypes
samples = vcf.columns[9:]
for s in samples:
    vcf[[f'{s}_1', f'{s}_2']] = vcf[s].str.split('|', expand=True)
vcf.drop(samples, inplace=True, axis=1)
vcf.columns = range(vcf.shape[1])

# read bed file
bed = pd.read_csv(bed_file, sep='\t', names=["chrom", "start", "end"])

dmax_list = np.zeros(bed.shape[0])

# iterate over all windows in bed
for i, row in bed.iterrows():
    print(f"Processing window {i+1} of {bed.shape[0]}", end='\r')
    wc = bed.iloc[i].iloc[0]
    ws = bed.iloc[i].iloc[1]
    we = bed.iloc[i].iloc[2]
    
    filtered_vcf = vcf[(vcf[0] == wc) & (vcf[1] >= ws) & (vcf[1] <= we)]
    hap_list = np.zeros(nhaps, dtype=object)
    for n in range(nhaps):
        hap_list[n] = ''.join(filtered_vcf.iloc[:, 9 + n].astype(str).values)
    
    d_list = np.zeros(int(len(hap_list) * (len(hap_list) - 1) / 2))
    for n in range(len(hap_list)):
        for m in range(n+1, len(hap_list)):
            d_list[int((n * len(hap_list) - n * (n + 1) / 2) + m - n - 1)] = ham_dist(hap_list[n], hap_list[m])

    dmax_list[i] = max(d_list)

bed['dmax'] = dmax_list
bed.to_csv(args.obed, sep='\t', index=False, header=False)
