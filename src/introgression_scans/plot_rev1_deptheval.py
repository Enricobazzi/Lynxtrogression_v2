import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

def get_files(idir, pop_pair, models):
    """
    Get the predictions files for the models to evaluate.
    """
    files = [f for f in os.listdir(idir) if pop_pair in f]
    files = [f for f in files if any(model in f for model in models)]
    files = [f for f in files if f.endswith('.predictions.csv')]
    return files

def get_preds_df_binary(idir, files, pthresh):
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

def get_preds_df_direction(idir, files, pthresh):
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
        if mig != 'none':
            df_list.append(df)
    
    preds_df = pd.concat([df for df in df_list if not df.empty])
    # apply my criteria to get prediction
    preds_df["ab_pred"] = np.where((preds_df["p_ab"] + preds_df["p_bi"] > pthresh), 1, 0)
    preds_df["ba_pred"] = np.where((preds_df["p_ba"] + preds_df["p_bi"] > pthresh), 1, 0)
    preds_df["bi_pred"] = np.where(preds_df["ab_pred"] + preds_df["ba_pred"] == 2, 1, 0)
    preds_df["none_pred"] = np.where((preds_df["ab_pred"] + preds_df["ba_pred"] + preds_df["bi_pred"]) == 0, 1, 0)
    preds_df["pred"] = np.where(preds_df["bi_pred"] == 1, 2, np.where(preds_df["ab_pred"] == 1, 0, np.where(preds_df["ba_pred"] == 1, 1, 3)))
    # remove predictions 3
    preds_df = preds_df[preds_df["pred"] != 3]
    # subsample to have 2000 rows of truth == 0 and 2000 rows of truth == 1 and 2000 rows of truth == 2
    preds_df_sub = pd.DataFrame()
    for truth in [0, 1, 2]:
        df_sub = preds_df[preds_df["truth"] == truth]
        if len(df_sub) > 2000:
            df_sub = df_sub.sample(n=2000, random_state=1)
        preds_df_sub = pd.concat([preds_df_sub, df_sub])
    preds_df = preds_df_sub
    preds_df = preds_df.reset_index(drop=True)
    return preds_df

def get_confusion_matrix(preds_df, indexes):
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
    cm = pd.DataFrame(cm, index=indexes, columns=indexes)
    cm.index.name = 'Simulated'
    cm.columns.name = 'Predicted'
    return cm, annot

def get_prec(cm, pred):
    tp = cm[pred].loc[pred]
    fp = cm[pred].sum() - tp
    return tp / (tp + fp)

def get_recall(cm, pred):
    tp = cm[pred].loc[pred]
    fn = cm.loc[pred].sum() - tp
    return tp / (tp + fn)

idir_base = 'data/introgression_scans/revision1'
pop_pairs = ['lpa-wel', 'lpa-sel']
folders = ['gterror_2X', 'gterror_3X', 'gterror_4X', 'gterror_5X', 'gterror_6X', 'gterror', 'evaluation_withM']
depths = ['2X', '3X', '4X', '5X', '6X', 'realistic', 'ideal']
pthresh = 0.9

for pop_pair in pop_pairs:
    if pop_pair == 'lpa-wel':
        models = ['12_9', '6_2', '20_7']
        pair = 'ILa-ELw'
    elif pop_pair == 'lpa-sel':
        models = ['12_6', '18_7', '18_10']
        pair = 'ILa-ELs'

    intro_p = []
    intro_r = []
    ab_p = []
    ab_r = []
    ba_p = []
    ba_r = []
    bi_p = []
    bi_r = []
    
    for folder, depth in zip(folders, depths):
        print(f'Processing depth: {depth}')
        idir = f'{idir_base}/{folder}'
        files = get_files(idir, pop_pair, models)
        preds_df = get_preds_df_binary(idir, files, pthresh=pthresh)
        cm, annot = get_confusion_matrix(preds_df, indexes=['intro', 'none'])
        prec_intro = get_prec(cm, 'intro')
        recall_intro = get_recall(cm, 'intro')
        preds_df = get_preds_df_direction(idir, files, pthresh=pthresh)
        cm, annot = get_confusion_matrix(preds_df, indexes=['ab', 'ba', 'bi'])
        prec_ab = get_prec(cm, 'ab')
        recall_ab = get_recall(cm, 'ab')
        prec_ba = get_prec(cm, 'ba')
        recall_ba = get_recall(cm, 'ba')
        prec_bi = get_prec(cm, 'bi')
        recall_bi = get_recall(cm, 'bi')
        intro_p.append(prec_intro)
        intro_r.append(recall_intro)
        ab_p.append(prec_ab)
        ab_r.append(recall_ab)
        ba_p.append(prec_ba)
        ba_r.append(recall_ba)
        bi_p.append(prec_bi)
        bi_r.append(recall_bi)

    # plot precision and recall for intro in one plot and 
    plt.figure(figsize=(6, 5))
    plt.plot(depths, intro_p, label='Precision', marker='o')
    plt.plot(depths, intro_r, label='Recall', marker='o')
    plt.ylim(0.7, 1)
    plt.xlabel('Evaluation dataset')
    plt.ylabel('Value')
    plt.title(f'{pair}\nIntrogression')
    plt.legend()
    plt.savefig(f'plots/revision1/{pair}_introgression_binary_prec_recall.png')
    
    # plot precision and recall for ab, ba, bi in two plots one for precision and one for recall
    fig, ax = plt.subplots(1, 2, figsize=(9, 5))
    ax[0].plot(depths, ab_p, label='ELtoIL', color='red', marker='o')
    ax[0].plot(depths, ba_p, label='ILtoEL', color='green', marker='o')
    ax[0].plot(depths, bi_p, label='BiDir', color='purple', marker='o')
    ax[0].set_ylim(0.4, 1)
    ax[0].set_xlabel('Evaluation dataset')
    ax[0].set_ylabel('Precision')
    ax[0].set_title(f'{pair}\ndirectionality precision')
    ax[0].legend()
    ax[1].plot(depths, ab_r, label='ELtoIL', color='red', marker='o')
    ax[1].plot(depths, ba_r, label='ILtoEL', color='green', marker='o')
    ax[1].plot(depths, bi_r, label='BiDir', color='purple', marker='o')
    ax[1].set_ylim(0.4, 1)
    ax[1].set_xlabel('Evaluation dataset')
    ax[1].set_ylabel('Recall')
    ax[1].set_title(f'{pair}\n directionality recall')
    ax[1].legend()
    plt.savefig(f'plots/revision1/{pair}_introgression_direction_multiclass_prec_recall.png')