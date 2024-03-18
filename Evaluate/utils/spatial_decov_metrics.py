"""
Adpoted from Li et al 2022
original file can be found at https://github.com/QuKunLab/SpatialBenchmarking/issues/6
"""
import os
import pandas as pd
import numpy as np
from collections.abc import Iterable
from sklearn.metrics import mean_squared_error
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, ttest_ind, mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns

def ssim(im1, im2, M=1):
    im1, im2 = im1 / im1.max(), im2 / im2.max()
    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1 = np.sqrt(((im1 - mu1) ** 2).mean())
    sigma2 = np.sqrt(((im2 - mu2) ** 2).mean())
    sigma12 = ((im1 - mu1) * (im2 - mu2)).mean()
    k1, k2, L = 0.01, 0.03, M
    C1 = (k1 * L) ** 2
    C2 = (k2 * L) ** 2
    C3 = C2 / 2
    l12 = (2 * mu1 * mu2 + C1) / (mu1**2 + mu2**2 + C1)
    c12 = (2 * sigma1 * sigma2 + C2) / (sigma1**2 + sigma2**2 + C2)
    s12 = (sigma12 + C3) / (sigma1 * sigma2 + C3)
    ssim = l12 * c12 * s12
    return ssim


def rmse(x1, x2):
    return mean_squared_error(x1, x2, squared=False)


def mae(x1, x2):
    return np.mean(np.abs(x1 - x2))


def compare_results(gd, result_list, metric="pcc", columns=None, axis=1):
    if metric == "pcc":
        func = pearsonr
        r_ind = 0
    if metric == "mae":
        func = mae
        r_ind = None
    if metric == "jsd":
        func = jensenshannon
        r_ind = None
    if metric == "rmse":
        func = rmse
        r_ind = None
    if metric == "ssim":
        func = ssim
        r_ind = None
    if isinstance(result_list, pd.DataFrame):
        c_list = []
        if axis == 1:
            print("axis: ", 1)
            for i, c in enumerate(gd.columns):
                r = func(gd.iloc[:, i].values, np.clip(result_list.iloc[:, i], 0, 1))
                if isinstance(result_list, Iterable):
                    if r_ind is not None:
                        r = r[r_ind]
                c_list.append(r)
        else:
            print("axis: ", 0)
            for i, c in enumerate(gd.index):
                r = func(gd.iloc[i, :].values, np.clip(result_list.iloc[i, :], 0, 1))
                if isinstance(result_list, Iterable):
                    if r_ind is not None:
                        r = r[r_ind]
                c_list.append(r)
        df = pd.DataFrame(c_list, index=gd.columns, columns=columns)
    else:
        df_list = []
        for res in result_list:
            c_list = []
            if axis == 1:
                for i, c in enumerate(gd.columns):
                    r = func(gd.iloc[:, i].values, np.clip(res.iloc[:, i], 0, 1))
                    if isinstance(res, Iterable):
                        if r_ind is not None:
                            r = r[r_ind]
                    c_list.append(r)
                df_tmp = pd.DataFrame(c_list, index=gd.columns)
            else:
                for i, c in enumerate(gd.index):
                    r = func(gd.iloc[i, :].values, np.clip(res.iloc[i, :], 0, 1))
                    if isinstance(res, Iterable):
                        if r_ind is not None:
                            r = r[r_ind]
                    c_list.append(r)
                df_tmp = pd.DataFrame(c_list, index=gd.index)
            df_list.append(df_tmp)
        df = pd.concat(df_list, axis=1)
        df.columns = columns
    return df

def make_score(gd_results, results_dict:dict, saving_path, scale_factor = 2):
    from pathlib import Path
    path = Path(saving_path)
    path.mkdir(parents=True, exist_ok=True)
    
    methods, results = list(results_dict.keys()), list(results_dict.values())
    
    Tools_score=[]
    score_col=[]
    list_up = list(range(1,len(methods)+1))
    list_down  = list(range(len(methods),0,-1))
    
    score_spots_pcc = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=0,
        metric='pcc'
    )
    score_spots_ssim = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=0,
        metric='ssim'
    )
    score_spots_rmse = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=0,
        metric='rmse'
    )
    score_spots_jsd = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=0,
        metric='jsd'
    )
    
    score_clusters_pcc = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=1,
        metric='pcc'
    )
    Tools_score.append(pd.Series(list_down, index=score_clusters_pcc.mean().sort_values(ascending=False).index))
    score_col.append('PCC')
    
    score_clusters_ssim = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=1,
        metric='ssim'
    )
    Tools_score.append(pd.Series(list_down, index=score_clusters_ssim.mean().sort_values(ascending=False).index))
    score_col.append('SSIM')
    
    score_clusters_rmse = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=1,
        metric='rmse'
    )
    Tools_score.append(pd.Series(list_up, index=score_clusters_rmse.mean().sort_values(ascending=False).index))
    score_col.append('RMSE')
    
    score_clusters_jsd = compare_results(
        gd_results,
        results,
        columns = methods,
        axis=1,
        metric='jsd'
    )
    Tools_score.append(pd.Series(list_up, index=score_clusters_jsd.mean().sort_values(ascending=False).index))
    score_col.append('JS')
    
    score=pd.concat([m for m in Tools_score],axis=1)
    score.columns = score_col
    score = score/len(methods)
    score['AS'] = score.mean(axis=1)
    score.sort_values('AS', ascending=False, inplace=True)
    score.to_csv(os.path.join(saving_path, 'AS_scores.csv'))
    score_clusters_pcc.to_csv(os.path.join(saving_path, 'PCC_scores.csv'))
    score_clusters_jsd.to_csv(os.path.join(saving_path, 'JSD_scores.csv'))
    score_clusters_rmse.to_csv(os.path.join(saving_path, 'RMSE_scores.csv'))
    score_clusters_ssim.to_csv(os.path.join(saving_path, 'SSIM_scores.csv'))

    #######
    # plot
    
    fig,axes = plt.subplots(ncols=4,nrows=2, figsize = (8 * scale_factor, 4 * scale_factor))
    order = list(score.index)
    colors = 'Set2' #["#9de846", '#F9EC31', "#BBA8D1","#D6DE23", "#7BD1F1"]
    sns.boxplot(data=score_spots_pcc,order=order,palette=colors,ax=axes[0,0],orient='h', showfliers = False,showmeans=True)
    axes[0,0].set_title('PCC of locations')
    sns.boxplot(data=score_spots_ssim,order=order,palette=colors,ax=axes[0,1],orient='h', showfliers = False,showmeans=True)
    axes[0,1].set_title('SSIM of locations')
    sns.boxplot(data=score_spots_rmse,order=order,palette=colors,ax=axes[0,2],orient='h', showfliers = False,showmeans=True)
    axes[0,2].set_title('RMSE of locations')
    sns.boxplot(data=score_spots_jsd,order=order,palette=colors,ax=axes[0,3],orient='h', showfliers = False,showmeans=True)
    axes[0,3].set_title('JSD of locations')
    
    sns.boxplot(data=score_clusters_pcc,order=order,palette=colors,ax=axes[1,0],orient='h', showfliers = False,showmeans=True)
    axes[1,0].set_title('PCC of clusters')
    sns.boxplot(data=score_clusters_ssim,order=order,palette=colors,ax=axes[1,1],orient='h', showfliers = False,showmeans=True)
    axes[1,1].set_title('SSIM of clusters')
    sns.boxplot(data=score_clusters_rmse,order=order,palette=colors,ax=axes[1,2],orient='h', showfliers = False,showmeans=True)
    axes[1,2].set_title('RMSE of clusters')
    sns.boxplot(data=score_clusters_jsd,order=order,palette=colors,ax=axes[1,3],orient='h', showfliers = False,showmeans=True)
    axes[1,3].set_title('JSD of clusters')
    
    for i in range(1,4):
        for j in [0,1]:
            axes[j,i].set_yticks([])
    
    # fig.tight_layout()
    fig.savefig(os.path.join(saving_path,'score_all_metrics.pdf'), bbox_inches='tight')
    plt.show()
    plt.close()