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

# def get_prec(df, pred):
#     return df[(df['truth'] == pred) & (df['pred'] == pred)].shape[0] / df[df['pred'] == pred].shape[0]

# def get_recall(df, pred):
#     return df[(df['truth'] == pred) & (df['pred'] == pred)].shape[0] / df[df['truth'] == pred].shape[0]

idir = 'data/introgression_scans/evaluation_withM'
pop_pairs = ['lpa-wel', 'lpa-sel']

for pop_pair in pop_pairs:
    if pop_pair == 'lpa-wel':
        models = ['12_9', '6_2', '20_7']
        pair = 'ILa-ELw'
    elif pop_pair == 'lpa-sel':
        models = ['12_6', '18_7', '18_10']
        pair = 'ILa-ELs'

    files = get_files(idir, pop_pair, models)
    pts = [0.75, 0.8, 0.85, 0.9, 0.95]
    
    intro_p = []
    intro_r = []
    ab_p = []
    ab_r = []
    ba_p = []
    ba_r = []
    bi_p = []
    bi_r = []
    
    for pthresh in pts:
        preds_df_direction = get_preds_df_direction(idir, files, pthresh)
        preds_df_binary = get_preds_df_binary(idir, files, pthresh)
        cm_binary, annot_binary = get_confusion_matrix(preds_df_binary, ['intro', 'none'])
        cm_direction, annot_direction = get_confusion_matrix(preds_df_direction, ['ab', 'ba', 'bi'])
        intro_p.append(get_prec(cm_binary, 'intro'))
        intro_r.append(get_recall(cm_binary, 'intro'))
        ab_p.append(get_prec(cm_direction, 'ab'))
        ab_r.append(get_recall(cm_direction, 'ab'))
        ba_p.append(get_prec(cm_direction, 'ba'))
        ba_r.append(get_recall(cm_direction, 'ba'))
        bi_p.append(get_prec(cm_direction, 'bi'))
        bi_r.append(get_recall(cm_direction, 'bi'))
    
    # plot precision and recall for intro in one plot
    plt.figure(figsize=(6, 5))
    
    plt.plot(pts, intro_p, label='Precision')
    plt.plot(pts, intro_r, label='Recall')
    plt.scatter(pts, intro_p, label='')
    plt.scatter(pts, intro_r, label='')
    plt.title(f'{pair}\nIntrogression')
    plt.xlabel('pthresh')
    plt.yticks([0.75, 0.8, 0.85, 0.9, 0.95])
    plt.xticks([0.75, 0.8, 0.85, 0.9, 0.95])
    plt.legend()
    plt.savefig(f'plots/introgression_scans/evaluation_withM/{pop_pair}_precision_recall.binary.pdf', format='pdf', bbox_inches='tight')
    
    # plot precision and recall for ab, ba, bi in two plots one for precision and one for recall
    fig, ax = plt.subplots(1, 2, figsize=(6, 5))
    ax[0].plot(pts, ab_p, label='ELtoIL', color='red')
    ax[0].plot(pts, ba_p, label='ILtoEL', color='green')
    ax[0].plot(pts, bi_p, label='BiDir', color='purple')
    ax[0].scatter(pts, ab_p, label='', color='red')
    ax[0].scatter(pts, ba_p, label='', color='green')
    ax[0].scatter(pts, bi_p, label='', color='purple')
    ax[0].set_title(f'{pair}\ndirectionality precision')
    ax[0].set_xlabel('pthresh')
    ax[0].set_yticks([0.75, 0.8, 0.85, 0.9, 0.95])
    ax[0].set_xticks([0.75, 0.8, 0.85, 0.9, 0.95])
    ax[1].plot(pts, ab_r, label='ELtoIL', color='red')
    ax[1].plot(pts, ba_r, label='ILtoEL', color='green')
    ax[1].plot(pts, bi_r, label='BiDir', color='purple')
    ax[1].scatter(pts, ab_r, label='', color='red')
    ax[1].scatter(pts, ba_r, label='', color='green')
    ax[1].scatter(pts, bi_r, label='', color='purple')
    ax[1].set_title(f'{pair}\ndirectionality recall')
    ax[1].set_xlabel('pthresh')
    ax[1].set_yticks([0.75, 0.8, 0.85, 0.9, 0.95])
    ax[1].set_xticks([0.75, 0.8, 0.85, 0.9, 0.95])
    plt.legend()
    plt.savefig(f'plots/introgression_scans/evaluation_withM/{pop_pair}_precision_recall.direction.pdf', format='pdf', bbox_inches='tight')

# for pop_pair in pop_pairs:
#     if pop_pair == 'lpa-wel':
#         models = ['12_9', '6_2', '20_7']
#     elif pop_pair == 'lpa-sel':
#         models = ['12_6', '18_7', '18_10'] 
#     # get table with following columns: ptresh, category, precision, recall
#     # category: ab, ba, bi, none (0, 1, 2, 3)
#     # pthesh: 0.75, 0.8, 0.85, 0.9, 0.95
#     # precision: TP / (TP + FP)
#     # recall: TP / (TP + FN)
#     pthreshs = [0.75, 0.8, 0.85, 0.9, 0.95]
#     results = []
#     for pthresh in pthreshs:
#         files = get_files(idir, pop_pair, models)
#         preds_df = get_preds_df(idir, files, pthresh)
#         for pred in range(4):
#             results.append([pthresh, pred, get_prec(preds_df, pred), get_recall(preds_df, pred)])
#     results = pd.DataFrame(results, columns=['pthresh', 'category', 'precision', 'recall'])
#     
#     # plot precision and fnr for each category across pthreshs
#     fig, ax = plt.subplots(1, 2, figsize=(8, 6))
#     for i, category in enumerate(['ab', 'ba', 'bi', 'none']):
#         df = results[results['category'] == i]
#         ax[0].plot(np.array(df['pthresh']), np.array(df['precision']))
#         ax[0].scatter(np.array(df['pthresh']), np.array(df['precision']), label=category)
#         ax[1].plot(np.array(df['pthresh']), np.array(df['recall']))
#         ax[1].scatter(np.array(df['pthresh']), np.array(df['recall']), label=category)
#     ax[0].set_title('Precision')
#     ax[1].set_title('Recall')
#     ax[0].set_xlabel('pthresh')
#     ax[1].set_xlabel('pthresh')
#     # place legend outside of plot
#     plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
#     plt.savefig(f'plots/introgression_scans/evaluation_withM/{pop_pair}_precision_recall.pdf', format='pdf', bbox_inches='tight')

