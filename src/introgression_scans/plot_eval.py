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
            df['truth'] = 0
        elif mig == 'ba':
            df['truth'] = 1
        elif mig == 'baab' or mig == 'abba':
            df['truth'] = 2
        elif mig == 'none':
            df['truth'] = 3
        df_list.append(df)
    preds_df = pd.concat([df for df in df_list if not df.empty])
    # apply my criteria to get prediction
    preds_df["ab_pred"] = np.where((preds_df["p_ab"] + preds_df["p_bi"] > pthresh), 1, 0)
    preds_df["ba_pred"] = np.where((preds_df["p_ba"] + preds_df["p_bi"] > pthresh), 1, 0)
    preds_df["bi_pred"] = np.where(preds_df["ab_pred"] + preds_df["ba_pred"] == 2, 1, 0)
    preds_df["none_pred"] = np.where((preds_df["ab_pred"] + preds_df["ba_pred"] + preds_df["bi_pred"]) == 0, 1, 0)
    preds_df["pred"] = np.where(preds_df["bi_pred"] == 1, 2, np.where(preds_df["ab_pred"] == 1, 0, np.where(preds_df["ba_pred"] == 1, 1, 3)))

    return preds_df

def get_confusion_matrix(preds_df):
    """
    Get the confusion matrix for the predictions. This is copied from introUnets.
    """
    y_true = preds_df["truth"]
    y_pred = preds_df["pred"]
    cm = confusion_matrix(y_true, y_pred)
    cm_sum = np.sum(cm, axis=1, keepdims=True)
    cm_perc = cm / cm_sum.astype(float) * 100
    annot = np.empty_like(cm).astype(str)
    nrows, ncols = cm.shape
    for i in range(nrows):
        for j in range(ncols):
            c = cm[i, j]
            p = cm_perc[i, j]
            if i == j:
                s = cm_sum[i]
                annot[i, j] = '%.1f%%\n%d/%d' % (p, c, s)
            elif c == 0:
                annot[i, j] = ''
            else:
                annot[i, j] = '%.1f%%\n%d' % (p, c)
    cm = pd.DataFrame(cm, index=['ab', 'ba', 'bi', 'none'], columns=['ab', 'ba', 'bi', 'none'])
    cm.index.name = 'Simulated'
    cm.columns.name = 'Predicted'
    return cm, annot

def plot_cm(cm, annot, oplot):
    """
    Plot the confusion matrix.
    """
    plt.figure(figsize=(10, 7))
    sns.heatmap(cm, annot=annot, fmt='', cmap='viridis', xticklabels=['ab', 'ba', 'bi', 'none'], yticklabels=['ab', 'ba', 'bi', 'none'])
    plt.savefig(oplot)

def main():
    args = parse_args()
    files = get_files(args.idir, args.pop_pair, args.models.split(','))
    preds_df = get_preds_df(args.idir, files, args.pthresh)
    cm, annot = get_confusion_matrix(preds_df)
    plot_cm(cm, annot, args.oplot)

if __name__ == "__main__":
    main()
