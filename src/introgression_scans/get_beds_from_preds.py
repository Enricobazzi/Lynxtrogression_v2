"""
Transform the predictions from apply_disc_to_npz.py into bed files of ab, ba, and bi introgressed windows.

Based on a probability threshold: 0.95 by default. The probability of no introgression must be < 0.01 (this is hard coded).

Needs a VCF file to get start and end positions of the windows.

Will get predictions from all csv files (extension = *.predictions.csv) in the input directory.

Window size and step sizes: 128 and 64 by default.

Will output a global predictions table and 3 bed files (ab, ba, bi) including all chromosomes.

Args:
    --folder: input directory with CSV files of predictions and where output table and beds are written.
    --ivcf: input VCF file.
    --pthresh: probability threshold.
    --wsize: window size.
    --step: step size.
"""

import os
import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type = str, required = True)
    parser.add_argument("--ivcf", type = str, required = True)
    parser.add_argument("--pthresh", type = float, default = 0.95)
    parser.add_argument("--wsize", type = int, default = 128)
    parser.add_argument("--step", type = int, default = 64)
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

def make_predictions(merged, pthresh):
    pred = merged.copy()
    pred["ab_pred"] = np.where((pred["p_ab"] + pred["p_bi"] > pthresh) & (pred["p_none"] < 0.01), 1, 0)
    pred["ba_pred"] = np.where((pred["p_ba"] + pred["p_bi"] > pthresh) & (pred["p_none"] < 0.01), 1, 0)
    pred["bi_pred"] = np.where(pred["ab_pred"] + pred["ba_pred"] == 2, 1, 0)
    return pred

def write_files(predictions, folder, pthresh):
    n = str(pthresh).replace(".", "")
    predictions.to_csv(os.path.join(folder, f"predictions_table_{n}.csv"), index = False)
    ab = predictions[predictions["ab_pred"] == 1]
    ab = ab[["chrom", "start", "end"]]
    ab.to_csv(os.path.join(folder, f"ab_introgressed_{n}.bed"), sep = "\t", header = False, index = False)
    ba = predictions[predictions["ba_pred"] == 1]
    ba = ba[["chrom", "start", "end"]]
    ba.to_csv(os.path.join(folder, f"ba_introgressed_{n}.bed"), sep = "\t", header = False, index = False)
    bi = predictions[predictions["bi_pred"] == 1]
    bi = bi[["chrom", "start", "end"]]
    bi.to_csv(os.path.join(folder, f"bi_introgressed_{n}.bed"), sep = "\t", header = False, index = False)

def main():
    args = parse_args()
    window_df = get_windows_df(args.ivcf, args.wsize, args.step)
    preds_df = get_preds_df(args.folder)
    merged_df = merge_dfs(window_df, preds_df)
    predictions = make_predictions(merged_df, args.pthresh)
    write_files(predictions, args.folder, args.pthresh)
    
if __name__ == "__main__":
    main()

