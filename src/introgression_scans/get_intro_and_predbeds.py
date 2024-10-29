import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

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
    """
    Get a DataFrame with the predictions from the introgression scans.
    """
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
    """
    Merge the windows and predictions DataFrames.
    """
    merged = pd.merge(wdf, pdf, on = ["chrom", "window"], how = "inner")
    merged["wlen"] = merged["end"] - merged["start"]
    merged["ab_prop"] = merged["p_ab"] * merged["wlen"]
    merged["ba_prop"] = merged["p_ba"] * merged["wlen"]
    merged["bi_prop"] = merged["p_bi"] * merged["wlen"]
    merged["none_prop"] = merged["p_none"] * merged["wlen"]
    return merged

def make_predictions(merged, pthresh):
    """
    Make the predictions based on the probabilities and the threshold.

    The predictions are made as follows:
    - If the sum of the probabilities of AB and BI is greater than the threshold, the window is called AB.
    - If the sum of the probabilities of BA and BI is greater than the threshold, the window is called BA.
    - If both AB and BA are called, the window is called BI.
    - If none of the above conditions are met, the window is called none.
    Additionally, windows adjacent to introgressed windows that have up to 0.05 less than the threshold are also called introgressed.
    """
    pred = merged.copy()
    # apply threshold
    pred["ab_pred"] = np.where((pred["p_ab"] + pred["p_bi"] > pthresh), 1, 0)
    pred["ba_pred"] = np.where((pred["p_ba"] + pred["p_bi"] > pthresh), 1, 0)
    pred["bi_pred"] = np.where(pred["ab_pred"] + pred["ba_pred"] == 2, 1, 0)
    # find windows adjacent to introgressed windows that have up to 0.1 less than the threshold and call them intro as well
    for i in range(1, len(pred) - 1):
        if pred.loc[i - 1, "ab_pred"] == 1 or pred.loc[i + 1, "ab_pred"] == 1:
            pred.loc[i, "ab_pred"] = np.where((pred.loc[i, "p_ab"] + pred.loc[i, "p_bi"] > pthresh - 0.2) & (pred.loc[i, "ab_pred"] == 0), 1, pred.loc[i, "ab_pred"])
        if pred.loc[i - 1, "ba_pred"] == 1 or pred.loc[i + 1, "ba_pred"] == 1:
            pred.loc[i, "ba_pred"] = np.where((pred.loc[i, "p_ba"] + pred.loc[i, "p_bi"] > pthresh - 0.2) & (pred.loc[i, "ba_pred"] == 0), 1, pred.loc[i, "ba_pred"])
    pred["ab_pred"] = pred["ab_pred"] * pred["wlen"]
    pred["ba_pred"] = pred["ba_pred"] * pred["wlen"]
    pred["bi_pred"] = pred["bi_pred"] * pred["wlen"]
    pred["none_pred"] = pred["wlen"] - pred["ab_pred"] - pred["ba_pred"] - pred["bi_pred"]
    pred["none_pred"] = [0 if x < 0 else x for x in pred["none_pred"]]
    return pred

# load VCFs and predictions
ivcf_wel = "data/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-wel.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf"
ipreds_wel = "data/introgression_scans/lpa-wel_predictions_withM"
ivcf_sel = "data/lynxtrogression_v2.autosomic_scaffolds.filter4.lpa-sel.ps.phased.merged.concat.fixed.afan.rd_fil.variant.vcf"
ipreds_sel = "data/introgression_scans/lpa-sel_predictions_withM"
wsize = 128
step = 64
wdf_wel = get_windows_df(ivcf_wel, wsize, step)
wdf_sel = get_windows_df(ivcf_sel, wsize, step)
pdf_wel = get_preds_df(ipreds_wel)
pdf_sel = get_preds_df(ipreds_sel)
merged_wel = merge_dfs(wdf_wel, pdf_wel)
merged_sel = merge_dfs(wdf_sel, pdf_sel)

# make predictions
wel_lpa = make_predictions(merged_wel, 0.9)
sel_lpa = make_predictions(merged_sel, 0.9)
wel_lpa.to_csv("data/introgression_scans/lpa-wel_predictions_withM.csv", index = False)
sel_lpa.to_csv("data/introgression_scans/lpa-sel_predictions_withM.csv", index = False)

# print amount of introgression - [NOTE] this is incorrect, as windows are overlapping
# in lpa from wel
# print(f'proportion of wel in lpa: {wel_lpa["ab_pred"].sum() / wel_lpa["wlen"].sum()}')
# # in lpa from sel
# print(f'proportion of sel in lpa: {sel_lpa["ab_pred"].sum() / sel_lpa["wlen"].sum()}')
# # in wel from lpa
# print(f'proportion of lpa in wel: {wel_lpa["ba_pred"].sum() / wel_lpa["wlen"].sum()}')
# # in sel from lpa
# print(f'proportion of lpa in sel: {sel_lpa["ba_pred"].sum() / sel_lpa["wlen"].sum()}')

# write bed files
# ab windows in wel_lpa
ab_wel = wel_lpa.copy()[wel_lpa["ab_pred"] > 0]
ab_wel = ab_wel[["chrom", "start", "end"]].reset_index(drop = True)
ab_wel.to_csv("data/introgression_scans/bed_files/wel_to_lpa_intro.bed", sep = "\t", header = False, index = False)
# ab windows in sel_lpa
ab_sel = sel_lpa.copy()[sel_lpa["ab_pred"] > 0]
ab_sel = ab_sel[["chrom", "start", "end"]].reset_index(drop = True)
ab_sel.to_csv("data/introgression_scans/bed_files/sel_to_lpa_intro.bed", sep = "\t", header = False, index = False)
# ba windows in wel_lpa
ba_wel = wel_lpa.copy()[wel_lpa["ba_pred"] > 0]
ba_wel = ba_wel[["chrom", "start", "end"]].reset_index(drop = True)
ba_wel.to_csv("data/introgression_scans/bed_files/lpa_to_wel_intro.bed", sep = "\t", header = False, index = False)
# ba windows in sel_lpa
ba_sel = sel_lpa.copy()[sel_lpa["ba_pred"] > 0]
ba_sel = ba_sel[["chrom", "start", "end"]].reset_index(drop = True)
ba_sel.to_csv("data/introgression_scans/bed_files/lpa_to_sel_intro.bed", sep = "\t", header = False, index = False)
