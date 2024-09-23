"""
Make a dataset of random windows of a given size and step size from a given VCF file.

Args:
--ivcf: input VCF file.
--wsize: window size.
--step: step size.
--nwin: number of windows to generate.
--ofolder: output folder for the BED files.
--nbeds: number of BED files to generate.
"""
import argparse
import pandas as pd
import random
import os

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ivcf", type = str, required = True)
    parser.add_argument("--wsize", type = int, required = True)
    parser.add_argument("--step", type = int, required = True)
    parser.add_argument("--nwin", type = int, required = True)
    parser.add_argument("--ofolder", type = str, required = True)
    parser.add_argument("--nbeds", type = int, default = 1)
    return parser.parse_args()

def get_windows_df(ivcf, wsize, step):
    """
    Get a DataFrame with the chromosome, start, and end positions of the windows (0-based).
    """
    vcf_dict = {}
    with open(ivcf, "r") as f:
        for line in f:
            if not line.startswith("#"):
                chrom, pos = line.strip().split("\t")[:2]
                vcf_dict.setdefault(chrom, []).append(int(pos))
    c1, c2, c3, c4 = [], [], [], []
    for chr in vcf_dict.keys():
        c = 0
        for i in range(0, len(vcf_dict[chr]), step):
            if i + wsize > len(vcf_dict[chr]):
                break
            start = vcf_dict[chr][i:i + wsize][0]
            end = vcf_dict[chr][i:i + wsize][-1]
            c1.append(chr)
            c2.append(start - 1)
            c3.append(end)
            c4.append(c)
            c += 1
    df = pd.DataFrame({"window": c4, "chrom": c1, "start": c2, "end": c3})
    return df

def get_random_wins(df, nwin):
    """
    Generate a DataFrame by randomly selecting nwin windows from the input DataFrame.
    """
    random_wins = random.sample(range(len(df)), nwin)
    # sort by window number
    random_wins.sort()
    return df.iloc[random_wins]

def write_bed(df, obed):
    """
    Write the DataFrame to a BED file.
    """
    df = df[["chrom", "start", "end"]]
    df.to_csv(obed, sep = "\t", header = False, index = False)

def main():
    args = parse_args()
    df = get_windows_df(args.ivcf, args.wsize, args.step)
    for i in range(args.nbeds):
        random_df = get_random_wins(df, args.nwin)
        obed = os.path.join(args.ofolder, f"random_windows.{i}.bed")
        print(f"Writing {obed}", end = "\r")
        write_bed(random_df, obed)
    print()

if __name__ == "__main__":
    main()
