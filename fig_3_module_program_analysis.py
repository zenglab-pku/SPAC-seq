import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import NMF
import pickle
from scipy.stats import pearsonr
from tqdm import tqdm
import consensusclustering as cc
from sklearn.cluster import AgglomerativeClustering

from A_utils import combine_guide_replicates

def run_nmf_and_consensus_clustering(
    fdata: ad.AnnData,
    n_top_genes: int = 10000,
    nmf_components: int = 100,
    nmf_max_iter: int = 2000,
    nmf_model_path: str = './nmf_model.hvg.100.pkl',
    nmf_matrix_path: str = './nmf_matrix.hvg.100.npy',
    corr_matrix_path: str = './nmf_corr_matrix.hvg.100.npy',
    cc_model_path: str = './cc_model.hvg.100.pkl',
    min_clusters: int = 4,
    max_clusters: int = 20,
    n_resamples: int = 100,
    resample_frac: float = 0.8,
    k_param: str = 'n_clusters'
):

    sc.pp.highly_variable_genes(fdata, n_top_genes=n_top_genes, flavor='cell_ranger', inplace=True)
    hvg_data = fdata[:, fdata.var['highly_variable']].copy()

    cnt_matrix = hvg_data.X.toarray()
    model = NMF(n_components=nmf_components, verbose=1, max_iter=nmf_max_iter)
    cnt_matrix_trans = model.fit_transform(cnt_matrix)
    pickle.dump(model, open(nmf_model_path, 'wb'))
    np.save(nmf_matrix_path, cnt_matrix_trans)

    num_components = cnt_matrix_trans.shape[1]
    corr_matrix = np.zeros((num_components, num_components))
    for i in tqdm(range(num_components)):
        for j in range(i, num_components):
            corr_matrix[i, j] = pearsonr(cnt_matrix_trans[:, i], cnt_matrix_trans[:, j])[0]
    for i in tqdm(range(num_components)):
        for j in range(i):
            corr_matrix[i, j] = corr_matrix[j, i]
    np.save(corr_matrix_path, corr_matrix)

    cc_model = cc.ConsensusClustering(
        AgglomerativeClustering(),
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        n_resamples=n_resamples,
        resample_frac=resample_frac,
        k_param=k_param
    )
    cc_model.fit(corr_matrix)
    cc_model.best_k()
    pickle.dump(cc_model, open(cc_model_path, 'wb'))

    return {
        'nmf_model': model,
        'nmf_matrix': cnt_matrix_trans,
        'corr_matrix': corr_matrix,
        'cc_model': cc_model
    }

def compute_normalized_guide_cluster_df(
    fdata: ad.AnnData,
    gdata: ad.AnnData,
    guides_of_interest: list[str] = None,
    cluster_key: str = 'cluster_cellcharter_givenk',
    samples: tuple = ('1-1', '1-2', '2-1', '2-2', '3-1', '3-2')
):
    if guides_of_interest is None:
        guides_of_interest = gdata.var_names.tolist()

    fdata.obs['coord_x'] = fdata.obsm['spatial'][:, 0] // 100 * 100
    fdata.obs['coord_y'] = fdata.obsm['spatial'][:, 1] // 100 * 100

    gdata.obs['coord_x'] = gdata.obsm['spatial'][:, 0] // 100 * 100
    gdata.obs['coord_y'] = gdata.obsm['spatial'][:, 1] // 100 * 100

    fdata.obs_names = [
        f"{marker}_{x}_{y}"
        for marker, x, y in zip(fdata.obs['marker'], fdata.obs['coord_x'], fdata.obs['coord_y'])
    ]
    gdata.obs_names = [
        f"{marker}_{x}_{y}"
        for marker, x, y in zip(gdata.obs['marker'], gdata.obs['coord_x'], gdata.obs['coord_y'])
    ]

    common_bins = np.intersect1d(fdata.obs_names, gdata.obs_names)
    fdata = fdata[common_bins].copy()
    gdata = gdata[common_bins].copy()

    gdata = combine_guide_replicates(gdata)
    gdata.obs['cluster'] = fdata.obs[cluster_key]

    df_list = []
    for sample in samples:
        mask_sample = gdata.obs.marker == sample

        df_mat = gdata[mask_sample, guides_of_interest].X
        index = gdata[mask_sample, guides_of_interest].obs_names
        df = pd.DataFrame(df_mat, index=index, columns=guides_of_interest)

        cluster_series = gdata[mask_sample, :].obs.cluster
        df = pd.concat([df, cluster_series], axis=1)
        df.cluster = df.cluster.astype(str)

        df = df.melt(
            id_vars=['cluster'],
            value_vars=guides_of_interest,
            var_name='guide',
            value_name='value'
        )
        df.value = df.value.astype(float)

        guide_sums = df.groupby('cluster')['value'].sum()
        df_agg = df.groupby(['cluster', 'guide'])['value'].sum().reset_index()
        df_agg['value'] = df_agg.apply(lambda x: x['value'] / guide_sums[x['cluster']], axis=1)
        df_agg['sample'] = sample
        df_list.append(df_agg)

    df_all = pd.concat(df_list)
    return df_all

def compute_guide_program_correlations(
    df_all,
    fdata,
    samples=['1-1', '1-2', '2-1', '2-2', '3-1', '3-2'],
    module_dict=None
):
    if module_dict is None:
        module_dict = {
            'Module 1': ['sgPomt1', 'sgCd44', 'sgMcoln1', 'sgS100a11'],
            'Module 2': ['sgDdit3', 'sgKcna3', 'sgCd52', 'sgL3mbtl3', 'sgAdrb2', 'sgSorl1', 'sgZhx2', 'sgWipf1', 'sgZscan12', 'sgTmem64', 'sgArntl', 'sgKlrd1'],
            'Module 3': ['sgNmb', 'sgPpia', 'sgScamp4', 'sgZc3h12a', 'sgFbxo7', 'sgCxcl16', 'sgIkbip', 'sgGpa33', 'sgPiezo1', 'sgSrgn', 'sgMark3', 'sgSlc39a8', 'sgErgic2', 'sgGata3', 'sgGlb1l2', 'sgAqp3', 'sgFlot1'],
        }
    plot_df_all = []
    for sample in samples:
        df_all_sample = df_all.loc[df_all['sample'] == sample, :].copy()

        nmf_scores = fdata.obs.groupby('cluster_cellcharter_givenk')[['Mask 0', 'Mask 1', 'Mask 2', 'Mask 3', 'Mask 4', 'Mask 5']].mean()

        nmf_scores_normalized = nmf_scores.apply(lambda x: (x - x.mean()) / np.sqrt(x.var()))

        results = []
        for cluster in df_all_sample.cluster.unique():
            cluster_data = df_all_sample[df_all_sample.cluster == cluster]
            
            cluster_nmf_scores = nmf_scores_normalized.loc[int(cluster)]
            
            for guide in cluster_data.guide.unique():
                guide_value = cluster_data[cluster_data.guide == guide].value.values[0]
                
                for i in range(6):
                    nmf_score = cluster_nmf_scores[f'Mask {i}']
                    correlation = guide_value * nmf_score      
                    results.append({
                        'cluster': cluster,
                        'guide': guide,
                        'program': i,
                        'correlation_score': correlation
                    })

        correlation_df = pd.DataFrame(results)
        plot_df = correlation_df.groupby(['program', 'guide'])['correlation_score'].mean().unstack()
        plot_df = plot_df.melt(ignore_index=False).reset_index(drop=False)
        plot_df['sample'] = sample

        plot_df_all.append(plot_df)

    plot_df_all = pd.concat(plot_df_all)
    
    plot_df_all['module'] = plot_df_all['guide'].apply(
        lambda x: next((module for module, guides in module_dict.items() if x in guides), None)
    )
    average_scores = plot_df_all.groupby(['module', 'program', 'sample'])['value'].mean().reset_index()

    sample = samples[-1]
    network_df = average_scores.loc[average_scores['sample'] == sample].copy()
    network_df.program = network_df.program.map(lambda x: 'Program ' + str(x + 1))

    return plot_df_all, average_scores, network_df