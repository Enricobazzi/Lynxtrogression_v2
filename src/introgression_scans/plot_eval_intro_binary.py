import os
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--idir", required = True, help = "Directory with predictions csv files (from apply_disc_to_npz.py).")
    parser.add_argument("--pop_pair", required = True, help = "Population pair as in the --pop-pair argument in apply_disc_to_npz.py.")
    parser.add_argument("--models", required = True, help = "Comma-separated list of model names to evaluate.")
    parser.add_argument("--pthresh", type = float, required=True, help = "Threshold for the predictions.")
    parser.add_argument("--oplot", required = True, help = "Output plot name.")
    return parser.parse_args()

def get_files(idir, pop_pair, models):
    """
    Get the predictions files for the models to evaluate.
    """
    files = [f for f in os.listdir(idir) if pop_pair in f]
    files = [f for f in files if any(model in f for model in models)]
    files = [f for f in files if f.endswith('.predictions.csv')]
    return files

def get_preds_df(idir, files, pthresh):
    """
    Get a DataFrame with the predictions from the introgression scans.
    """
    df_list = []
    for file in files:
        mig = file.split('.')[0].split('_')[-1]
        df = pd.read_csv(f'{idir}/{file}')
        if mig == 'ab':
            # select 333 random rows
            df = df.sample(n=333, random_state=1)
            df['truth'] = 0
        elif mig == 'ba':
            # select 333 random rows
            df = df.sample(n=333, random_state=1)
            df['truth'] = 0
        elif mig == 'baab' or mig == 'abba':
            # select 167 random rows
            df = df.sample(n=167, random_state=1)
            df['truth'] = 0
        elif mig == 'none':
            df['truth'] = 1
        df_list.append(df)
    
    preds_df = pd.concat([df for df in df_list if not df.empty])
    # apply my criteria to get prediction
    preds_df["ab_pred"] = np.where((preds_df["p_ab"] + preds_df["p_bi"] > pthresh), 1, 0)
    preds_df["ba_pred"] = np.where((preds_df["p_ba"] + preds_df["p_bi"] > pthresh), 1, 0)
    preds_df["bi_pred"] = np.where(preds_df["ab_pred"] + preds_df["ba_pred"] == 2, 1, 0)
    preds_df["none_pred"] = np.where((preds_df["ab_pred"] + preds_df["ba_pred"] + preds_df["bi_pred"]) == 0, 1, 0)
    preds_df["pred"] = np.where(preds_df["bi_pred"] == 1, 0, np.where(preds_df["ab_pred"] == 1, 0, np.where(preds_df["ba_pred"] == 1, 0, 1)))
    return preds_df

def get_confusion_matrix(preds_df, normalize="precision"):
    """
    Args:
        preds_df: DataFrame with 'truth' and 'pred' columns
        normalize: "recall" or "precision" (default)
    Returns:
        cm: confusion matrix as DataFrame
        annot: annotations to plot
    """
    y_true = preds_df["truth"]
    y_pred = preds_df["pred"]
    cm = confusion_matrix(y_true, y_pred)
    
    if normalize == "recall":
        cm_sum = np.sum(cm, axis=1, keepdims=True)  # normalize across rows
    elif normalize == "precision":
        cm_sum = np.sum(cm, axis=0, keepdims=True)  # normalize across columns
    else:
        raise ValueError("normalize must be 'recall' or 'precision'")
    
    cm_perc = cm / cm_sum.astype(float) * 100
    
    annot = np.empty_like(cm).astype(str)
    nrows, ncols = cm.shape
    for i in range(nrows):
        for j in range(ncols):
            c = cm[i, j]
            p = cm_perc[i, j]
            if i == j:
                s = cm_sum[0, j] if normalize == "precision" else cm_sum[i, 0]
                annot[i, j] = '%.1f%%\n%d/%d' % (p, c, s)
            elif c == 0:
                annot[i, j] = ''
            else:
                annot[i, j] = '%.1f%%\n%d' % (p, c)

    cm = pd.DataFrame(cm, index=['intro', 'none'], columns=['intro', 'none'])
    cm.index.name = 'Simulated'
    cm.columns.name = 'Predicted'
    return cm, annot

def plot_cm(cm, annot, oplot):
    plt.figure(figsize=(7, 4.9))
    # Convert counts to percentages for coloring
    cm_percent = cm.div(cm.sum(axis=0), axis=1) * 100  # axis=0 for precision coloring
    # If you want to color by recall instead, you would use axis=1
    sns.heatmap(cm_percent, annot=annot, fmt='', cmap='viridis', 
                xticklabels=['intro', 'none'], yticklabels=['intro', 'none'],
                vmin=0, vmax=100, cbar_kws={'label': 'Precision (%)'})
    
    plt.savefig(oplot)
    plt.close()

def print_summary(cm):
    tp = cm['intro'].iloc[0]
    tn = cm['none'].iloc[1]
    fp = cm['none'].iloc[0]
    fn = cm['intro'].iloc[1]
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    print(f'Precision: {precision:.3f}')
    print(f'Recall: {recall:.3f}')

def main():
    args = parse_args()
    files = get_files(args.idir, args.pop_pair, args.models.split(','))
    preds_df = get_preds_df(args.idir, files, args.pthresh)
    cm, annot = get_confusion_matrix(preds_df)
    plot_cm(cm, annot, args.oplot)
    print(f"{args.pop_pair} precision and recall:")
    print_summary(cm)

if __name__ == "__main__":
    main()
