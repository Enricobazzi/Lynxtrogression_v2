import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

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

def get_prec(df, pred):
    return df[(df['truth'] == pred) & (df['pred'] == pred)].shape[0] / df[df['pred'] == pred].shape[0]

def get_recall(df, pred):
    return df[(df['truth'] == pred) & (df['pred'] == pred)].shape[0] / df[df['truth'] == pred].shape[0]

idir = 'data/introgression_scans/evaluation_withM'
pop_pairs = ['lpa-wel', 'lpa-sel']

for pop_pair in pop_pairs:
    if pop_pair == 'lpa-wel':
        models = ['12_9', '6_2', '20_7']
    elif pop_pair == 'lpa-sel':
        models = ['12_6', '18_7', '18_10'] 
    # get table with following columns: ptresh, category, precision, recall
    # category: ab, ba, bi, none (0, 1, 2, 3)
    # pthesh: 0.75, 0.8, 0.85, 0.9, 0.95
    # precision: TP / (TP + FP)
    # recall: TP / (TP + FN)
    pthreshs = [0.75, 0.8, 0.85, 0.9, 0.95]
    results = []
    for pthresh in pthreshs:
        files = get_files(idir, pop_pair, models)
        preds_df = get_preds_df(idir, files, pthresh)
        for pred in range(4):
            results.append([pthresh, pred, get_prec(preds_df, pred), get_recall(preds_df, pred)])
    results = pd.DataFrame(results, columns=['pthresh', 'category', 'precision', 'recall'])
    
    # plot precision and fnr for each category across pthreshs
    fig, ax = plt.subplots(1, 2, figsize=(8, 6))
    for i, category in enumerate(['ab', 'ba', 'bi', 'none']):
        df = results[results['category'] == i]
        ax[0].plot(np.array(df['pthresh']), np.array(df['precision']))
        ax[0].scatter(np.array(df['pthresh']), np.array(df['precision']), label=category)
        ax[1].plot(np.array(df['pthresh']), np.array(df['recall']))
        ax[1].scatter(np.array(df['pthresh']), np.array(df['recall']), label=category)
    ax[0].set_title('Precision')
    ax[1].set_title('Recall')
    ax[0].set_xlabel('pthresh')
    ax[1].set_xlabel('pthresh')
    # place legend outside of plot
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.savefig(f'plots/introgression_scans/evaluation_withM/{pop_pair}_precision_recall.pdf', format='pdf', bbox_inches='tight')