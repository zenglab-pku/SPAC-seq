import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from tqdm import tqdm

def gene_mouse(fdata, str_):
    """
    Usage:
        Return mouse version of genes, check for absent genes
    """
    str_list = [t[0] + t[1:].lower() for t in str_]
    for str_ in str_list:
        if str_ not in fdata.var_names: print("Warning: " + str_ + " not in data variables!")
    return str_list

def remove_mito_ribo_hk_lnc_genes(adata, housekeeping_list="He2020Nature_mouseHK.txt"):
    
    return_data = adata.copy()
    return_data.var["mt"] = return_data.var_names.str.startswith("Mt")
    return_data.var["mt-"] = return_data.var_names.str.startswith("mt-")
    return_data.var["gm"] = return_data.var_names.str.startswith("Gm")
    return_data.var["Rb"] = return_data.var_names.str.startswith("Rp")
    return_data.var["rik"] = [True if "Rik" in str else False for str in return_data.var_names]
    return_data = return_data[:, ~return_data.var["mt"]].copy()
    return_data = return_data[:, ~return_data.var["mt-"]]
    return_data = return_data[:, ~return_data.var["Rb"]]
    return_data = return_data[:, ~return_data.var["gm"]]
    return_data = return_data[:, ~return_data.var["rik"]]

    with open(housekeeping_list, 'r') as f:
        for line in f:
            hk_genes = line.split('\t')
            break
    return_data = return_data[:, [gene for gene in return_data.var_names if gene not in hk_genes]]
    return return_data

def load_data_list_from_h5(data_list, marker_list=[], marker_name="marker"):
    """
    Usage:
        load file names <data_list> to anndata, concatenate them to one anndata
    Return:
        a object names unique anndata
    """
    anndata_list = []
    if len(marker_list) == 0:
        for idx, data in enumerate(data_list):
            anndata_ = sc.read_h5ad(data)
            anndata_list.append(anndata_)
    else:
        for idx, data in enumerate(data_list):
            anndata_ = sc.read_h5ad(data)
            anndata_.obs[marker_name] = marker_list[idx]
            anndata_list.append(anndata_)
    result_data = ad.concat(anndata_list)
    result_data.obs_names_make_unique()
    return result_data

def combine_guide_replicates(gdata):
    """
    Usage:
        combine guide replicates in <gdata> to a single gene name
    Return:
        single gene name anndata
    """
    sgs = gdata.var_names.str.split('_', n=1).str[0]
    sgs_grouped = pd.DataFrame(gdata.X.toarray(), columns=gdata.var_names)
    sgs_grouped = sgs_grouped.groupby(sgs, axis=1).sum()

    cgdata = ad.AnnData(sgs_grouped, obs=gdata.obs, var=pd.DataFrame(index=sgs_grouped.columns))
    cgdata.obsm['spatial'] = gdata.obsm['spatial']
    return cgdata

def calculate_deg(fdata, p_value_cuttoff=0.01, n_genes=200, ):
    """
    Usage:
        caculate deg from anndata <fdata>, sorts by score
        requires fdata to have perform rank_genes_groups() analysis
    Return:
        tuple (DEG list, A dataframe containing top degs)
    """
    result = fdata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    top_genes_df = pd.DataFrame()
    for group in groups:
        genes = result['names'][group]
        log2fc = result['logfoldchanges'][group]
        padj = result['pvals_adj'][group]
        score = result['scores'][group]
        data = pd.DataFrame({
            'Gene': genes,
            'Log2FoldChange': log2fc,
            'padj': padj,
            'score': score,
            'Cluster': group
        })
        filtered_data = data[data['padj'] < p_value_cuttoff]
        sorted_data = filtered_data.sort_values(by='score', ascending=False)
        top_genes = sorted_data.head(n_genes)
        top_genes_df = pd.concat([top_genes_df, top_genes], ignore_index=True)

    deg = top_genes_df["Gene"].unique().tolist()
    return deg, top_genes_df

def align_rna_guide_data(rnadata, guidedata, marker_name="marker"):
    """
    Usage:
        align <rnadata> and <guidedata> together, to ensure their <marker_name> field + spatial coordinate match
        note that <rnadata> is used as reference, any bin not in <guidedata> or excessive from <rnadata> is add or removed, respectively
        also sorts the two anndata by adding the 'obs': "cov"
    Return:
        a tuple, the sorted rna anndata and the cleaned sorted guide anndata
    """

    rnadata.obs["cov"] = [str(marker) + '_' + str(array[0]) + "-" + str(array[1]) for marker, array in zip(rnadata.obs[marker_name], rnadata.obsm["spatial"])]
    guidedata.obs["cov"] = [str(marker) + '_' + str(array[0]) + "-" + str(array[1]) for marker, array in zip(guidedata.obs[marker_name], guidedata.obsm["spatial"])]

    guidedata.obs_names = guidedata.obs['cov']
    rnadata.obs_names = rnadata.obs['cov']
    sorted_obs_names = rnadata.obs_names.sort_values()

    common_cov = np.intersect1d(rnadata.obs['cov'], guidedata.obs['cov'])
    guidedata_filtered = guidedata[common_cov].copy()
    guidedata_filtered = guidedata_filtered[sorted_obs_names].copy()
    rnadata_reordered = rnadata[sorted_obs_names].copy()

    missing_cov = rnadata.obs['cov'][~rnadata.obs['cov'].isin(common_cov)]
    if len(missing_cov) > 0:
        missing_obs = pd.DataFrame({marker_name: [cov.split('_')[0] for cov in missing_cov]}, index=missing_cov)
        missing_obsm_spatial = np.array([[float(cov.split('_')[1].split('-')[0]), float(cov.split('_')[1].split('-')[1])] for cov in missing_cov])
        missing_X = np.zeros((len(missing_cov), guidedata.shape[1]))

        missing_guidedata = ad.AnnData(X=missing_X, obs=missing_obs, obsm={'spatial': missing_obsm_spatial})
        guidedata_filtered = guidedata_filtered.concatenate(missing_guidedata)
    return rnadata_reordered, guidedata_filtered

def beta_regression(rnadata, guidedata, expr_term=[], expr_layer=None, \
    alpha_list=[0.001, 0.005, 0.01, 0.05, 0.1], max_iter=10000, l1_ratio=0.5):

    """
    Usage:
        beta regression
    Returns:
        dataframe of beta vector
    """

    from sklearn.linear_model import ElasticNet
    from sklearn.model_selection import GridSearchCV
    if len(rnadata.obs_names) != len(guidedata.obs_names):
        print("Error!, data bin count doesn't match!")
        return None
    if type(guidedata.X) == "scipy.sparse._csr.csr_matrix":
        guide_matrix = guidedata.X.toarray()
    else:
        guide_matrix = guidedata.X
    
    if all([obj in rnadata.obs for obj in expr_term]):
        print("Extracting bin metric")
        score_df = rnadata.obs[expr_term]
        try:
            score_matrix = np.array(score_df.astype(int))
        except:
            print("Error converting obs to np.float matrix!")
            return None
    elif all([obj in rnadata.var_names for obj in expr_term]):
        print("Extracting expression metric")
        if not expr_layer:
            if type(rnadata.X) == "scipy.sparse._csr.csr_matrix":
                score_matrix = rnadata[:, expr_term].X.toarray().T
            else:
                score_matrix = rnadata[:, expr_term].X.T
        else:
            if type(rnadata.layer[expr_layer]) == "scipy.sparse._csr.csr_matrix":
                score_matrix = rnadata[:, expr_term].layer[expr_layer].toarray().T
            else:
                score_matrix = rnadata[:, expr_term].layer[expr_layer].T
    
    corr_matrix = np.zeros((guide_matrix.shape[1], score_matrix.shape[1]))
    print("Regressing with " + str(guide_matrix.shape[0]) + " samples, " + str(guide_matrix.shape[1]) +\
        " guides, " + str(score_matrix.shape[1]) + " objectives.")
    param_grid = { 'alpha': alpha_list }
    range_ = range(guide_matrix.shape[1])
    for i in tqdm(range_):
        guide_vector = guide_matrix[:, i]
        regr=ElasticNet(max_iter=max_iter, l1_ratio=l1_ratio)
        grid_regr = GridSearchCV(regr, param_grid, scoring='neg_mean_squared_error', cv=5, n_jobs=-1, verbose=0)
        grid_regr.fit(X=score_matrix, y=guide_vector)
        regr = grid_regr.best_estimator_
        regr.fit(X=score_matrix, y=guide_vector)
        corr_matrix[i] = regr.coef_
    corr_df = pd.DataFrame(corr_matrix, index=guidedata.var_names, columns=expr_term)
    print("Regression done!")
    return corr_df

def calculate_correlation(rnadata, guidedata, expr_term=[], expr_layer=None, corr_type="spearmanr", control=False):
    """
    Usage:
        calculates the correlation with rndata <expr_term> and guidedata guidecount
        the <corr_type> function specifies the correlation method, can be "spearmanr" and "pearsonr"
    Return:
        the dataframe containing each rnadata <expr_term> and guide
    """

    from scipy.stats import pearsonr
    from scipy.stats import spearmanr
    if len(rnadata.obs_names) != len(guidedata.obs_names):
        print("Error!, data bin count doesn't match!")
        return None
    if type(guidedata.X) == "scipy.sparse._csr.csr_matrix":
        guide_matrix = guidedata.X.toarray()
    else:
        guide_matrix = guidedata.X
    if control:
        ntc_index = np.where(guidedata.var_names == "sgNTC")[0][0]
        ntc_vector = guide_matrix[:, ntc_index].reshape(-1, 1)
    
    if all([obj in rnadata.obs for obj in expr_term]):
        print("Extracting bin metric")
        score_df = rnadata.obs[expr_term]
        try:
            score_matrix = np.array(score_df.astype(int))
        except:
            print("Error converting obs to np.float matrix!")
            return None
    elif all([obj in rnadata.var_names for obj in expr_term]):
        print("Extracting expression metric")
        if not expr_layer:
            if type(rnadata.X) == "scipy.sparse._csr.csr_matrix":
                score_matrix = rnadata[:, expr_term].X.toarray().T
            else:
                score_matrix = rnadata[:, expr_term].X.T
        else:
            if type(rnadata.layer[expr_layer]) == "scipy.sparse._csr.csr_matrix":
                score_matrix = rnadata[:, expr_term].layer[expr_layer].toarray().T
            else:
                score_matrix = rnadata[:, expr_term].layer[expr_layer].T
    
    corr_matrix = np.zeros((guide_matrix.shape[1], score_matrix.shape[1]))
    print("Correlation with " + str(guide_matrix.shape[0]) + " samples, " + str(guide_matrix.shape[1]) +\
        " guides, " + str(score_matrix.shape[1]) + " objectives.")
    range_ = range(guide_matrix.shape[1])
    if control:
        if corr_type == "spearmanr":
            for i in tqdm(range_):
                if i == ntc_index: continue
                for j in range(score_matrix.shape[1]):
                    guide_vector = guide_matrix[:, i] - ntc_vector.T[0]
                    score_vector = score_matrix[:, j]
                    corr_matrix[i, j] = spearmanr(guide_vector, score_vector)[0]
        elif corr_type == "pearsonr":
            for i in tqdm(range_):
                if i == ntc_index: continue
                for j in range(score_matrix.shape[1]):
                    guide_vector = guide_matrix[:, i] - ntc_vector.T[0]
                    score_vector = score_matrix[:, j]
                    corr_matrix[i, j] = pearsonr(guide_vector, score_vector)[0]
        else:
            print("Error, cannot find correlation type " + corr_type)
            return None
    else:
        if corr_type == "spearmanr":
            for i in tqdm(range_):
                for j in range(score_matrix.shape[1]):
                    guide_vector = guide_matrix[:, i]
                    score_vector = score_matrix[:, j]
                    corr_matrix[i, j] = spearmanr(guide_vector, score_vector)[0]
        elif corr_type == "pearsonr":
            for i in tqdm(range_):
                for j in range(score_matrix.shape[1]):
                    guide_vector = guide_matrix[:, i]
                    score_vector = score_matrix[:, j]
                    corr_matrix[i, j] = pearsonr(guide_vector, score_vector)[0]
        else:
            print("Error, cannot find correlation type " + corr_type)
            return None
    corr_df = pd.DataFrame(corr_matrix, index=guidedata.var_names, columns=expr_term)
    if control:
        corr_df = corr_df.loc[corr_df.index != "sgNTC"]
    print("Regression done!")
    return corr_df

def extract_small_bin(bin20_data, bin100_masked_data, replace_obs=None, bin_size=100):
    """
    Usage:
        from big bin data extract corresponding spatial coordinate small bin data
        note that the data must be with the same chip, without hetero batches
    Return:
        extracted small bin size data
    """
    try:
        common_spatial_cov = bin20_data.obs['spatial_cov'].isin(bin100_masked_data.obs['spatial_cov']).values
    except:
        bin20_data.obs["spatial_cov"] = ['_'.join(map(str, [array[0] // bin_size * bin_size, array[1] // bin_size * bin_size])) for array in bin20_data.obsm["spatial"]]
        bin100_masked_data.obs["spatial_cov"] = ['_'.join(map(str, [array[0] // bin_size * bin_size, array[1] // bin_size * bin_size])) for array in bin100_masked_data.obsm["spatial"]]
        common_spatial_cov = bin20_data.obs['spatial_cov'].isin(bin100_masked_data.obs['spatial_cov']).values

    if replace_obs:
        for obj in replace_obs:
            spatial_to_obs = dict(zip(bin100_masked_data.obs["spatial_cov"], bin100_masked_data.obs[obj]))
            bin20_data.obs[obj] = bin20_data.obs["spatial_cov"].map(spatial_to_obs)
    return bin20_data[common_spatial_cov].copy()

def extract_big_bin(bin100_data, bin20_masked_data, bin_size=100):
    """
    Usage:
        from small bin data extract corresponding spatial coordinate big bin data
        note that the data must be with the same chip, without hetero batches
    Return:
        extracted big bin size data
    """
    try:
        bin20_mask = bin20_masked_data.obs["spatial_cov"].unique().tolist()
        mask = [array in bin20_mask for array in bin100_data.obs["spmatial_cov"]]
    except:
        bin100_data.obs["spatial_cov"] = ['_'.join(map(str, [array[0] // bin_size * bin_size, array[1] // bin_size * bin_size])) for array in bin100_data.obsm["spatial"]]
        bin20_masked_data.obs["spatial_cov"] = ['_'.join(map(str, [array[0] // bin_size * bin_size, array[1] // bin_size * bin_size])) for array in bin20_masked_data.obsm["spatial"]]
        bin20_mask = bin20_masked_data.obs["spatial_cov"].unique().tolist()
        mask = [array in bin20_mask for array in bin100_data.obs["spatial_cov"]]
    return bin100_data[mask].copy()

def find_near(bin100, bin100_data, border=1, bin_size=100):
    """
    Usage:
        extract bin with <border> number of bin(s) expanded from <bin100> in <bin100_data>
        note that the data must be with the same chip, without hetero batches
    Return:
        extracted neighbor anndata
    """
    cmin = ((bin100.obsm["spatial"] // bin_size - border) * bin_size)[0]
    cmax = ((bin100.obsm["spatial"] // bin_size + border + 1) * bin_size)[0]
    mask = [cmin[0] <= array[0] <= cmax[0] and cmin[1] <= array[1] <= cmax[1] for array in bin100_data.obsm["spatial"]]
    return bin100_data[mask]

def extract_region(fdata, x_region, y_region):
    df = pd.DataFrame(fdata.obsm["spatial"])
    mask = (df[0] > x_region[0]) & (df[0] < x_region[1]) & (df[1] > y_region[0]) & (df[1] < y_region[1])
    return fdata[mask].copy()

def plot_score(fdata, axs, marker_field="marker", field="score", size=2, cmap="viridis", vmax=None, vmin=None):
    """
    Usage:
        plot score with squidpy field <score>
        specify <marker_field> to draw different chips at the same time
    Return:
        various graph with <marker_field> as batch
    """
    if axs.shape[0] * axs.shape[1] != len(np.unique(fdata.obs[marker_field])):
        print("Error, chip size does not match ax size")
    for idx, marker in enumerate(np.unique(fdata.obs[marker_field])):
        ax = axs[idx // 2, idx % 2]
        sq.pl.spatial_scatter(fdata[fdata.obs[marker_field] == marker], color=field, size=size, shape=None, library_id="spatial", cmap=cmap, vmax=vmax, vmin=vmin, ax=ax)

def plot_contour(adata, score_name, ax, cmap="viridis", levels=8, fill=True):
    """
    Usage:
        plot contour map at pyplot <ax> axes with <score_name> in <adata>
    Return:
        a contour map
    """
    from scipy.interpolate import griddata

    apo_bin = adata.obs[score_name]
    hypo_coord = adata.obsm['spatial']
    coord_df = pd.DataFrame({
        'x_coord': hypo_coord[:, 0],
        'y_coord': hypo_coord[:, 1],
        'value': apo_bin
    })
    coord_df_piv = pd.pivot_table(coord_df, values="value", index="y_coord", columns="x_coord", aggfunc=np.mean)
    z = coord_df_piv.values
    y = coord_df_piv.index
    x = coord_df_piv.columns
    color_map = plt.get_cmap(cmap)
    print(z)

    new_cmap = mcolors.ListedColormap(color_map(np.linspace(0, 1, 100)))
    new_cmap.colors[:, -1] = np.linspace(0, 0.8, new_cmap.N)

    if not fill:
        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        xi, yi = np.meshgrid(xi, yi)

        zi = griddata((x.tolist(), y.tolist()), z.flatten(), (xi, yi), method='linear')
        #print(x.tolist(), y.tolist(), z, zi)
        x = xi
        y = yi
        z = zi
    
    cs = ax.contourf(x, y, z, levels=levels, cmap=new_cmap)
    plt.colorbar(cs, ax=ax, shrink=0.7)

def plot_dot(adata, gene_name, ax, color="Orange", square=True, alpha=0.2):
    """
    Usage:
        plot dot map at pyplot <ax> axes with <gene_name> in <adata>
        size indicate the <gene_name> count
    Return:
        a dot map
    """
    x = adata.obsm["spatial"][:, 0]
    y = adata.obsm["spatial"][:, 1]
    a = adata[:, gene_name].X.toarray().flatten()
    if square is True:
        ax.scatter(x, y, s=np.square(a), alpha=alpha, color=color, label=gene_name)
    else:
        ax.scatter(x, y, s=a, alpha=alpha, color=color, label=gene_name)