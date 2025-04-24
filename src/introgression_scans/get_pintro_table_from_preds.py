"""
Transform the predictions from apply_disc_to_npz.py into bed files of lpa-wel_pintro and lpa-sel_pintro windows table files.

The output table files will have the following format:
    window             chrom  start    end      p_ab      p_ba      p_bi    p_none

The positions are 0-based and the windows are numbered from 0 to n-1, where n is the number of windows.

Needs a VCF file to get start and end positions of the windows.

Window size and step sizes: 128 and 64 by default.
"""
import pandas as pd
import numpy as np
import os

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

def get_preds_df(ipreds):
    df_exists = None
    for file in os.listdir(ipreds):
        if file.endswith(".predictions.csv"):
            preds = pd.read_csv(os.path.join(ipreds, file))
            chr = '.'.join(file.split(".")[0:2])
            preds["chrom"] = chr
            if not df_exists:
                df = preds.copy()
                df_exists = True
            else:
                df = pd.concat([df, preds])
                df = df.reset_index(drop = True)
    return df

def merge_dfs(wdf, pdf):
    merged = pd.merge(wdf, pdf, on = ["chrom", "window"], how = "inner")
    return merged

wsize = 128
step = 64

for pop_pair in ['lpa-wel', 'lpa-sel']:
    folder = f'data/introgression_scans/{pop_pair}_predictions_withM'
    ivcf = f'data/lynxtrogression_v2.autosomic_scaffolds.filter4.{pop_pair}.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf'
    otable = f'{folder}/predictions_table.tsv'
    
    window_df = get_windows_df(ivcf, wsize, step)
    preds_df = get_preds_df(folder)
    merged_df = merge_dfs(window_df, preds_df)
    
    merged_df.to_csv(otable, sep = "\t", index = False)