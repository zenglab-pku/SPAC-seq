import anndata as ad
import scanpy as sc
import squidpy as sq
import cellcharter as cc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from anndata import AnnData
from typing import Optional
import sklearn.neighbors as sk_neighbors
import scipy.sparse as sp_sparse
import subprocess
from scipy.cluster.hierarchy import linkage
import scvi

from fast_utils import(
    compute_sseq_params_o3,
    sseq_differential_expression_o3
)

LOUVAIN_BIN_PATH = "louvain"
CONVERT_BIN_PATH = "convert"
LOUVAIN_PATH = "louvain_tmp"
LOUVAIN_DEFAULT_SEED = 0

def hdc(
    adata: AnnData,
    pca_key: str = "X_pca",
    leaf_size: int = 40,
    num_neighbors: int = 1,
    neighbor_a: int = -230,
    neighbor_b: int = 120,
    random_seed: int = LOUVAIN_DEFAULT_SEED,
    result_key: str = "hdc"
) -> Optional[AnnData]:

    expr_matrix = adata.X.T
    pca_matrix = adata.obsm[pca_key]
    num_cells = len(adata.obs_names)

    use_neighbors = int(max(num_neighbors, np.round(neighbor_a + neighbor_b * np.log10(num_cells))))
    n_neighbors = max(1, min(use_neighbors, num_cells - 1))

    balltree = sk_neighbors.BallTree(pca_matrix, leaf_size=leaf_size)
    nn_dist, nn_idx = balltree.query(pca_matrix, k=n_neighbors + 1)
    nn_idx = nn_idx[:, 1:]
    nn_dist = nn_dist[:, 1:]

    row_indices = np.repeat(np.arange(nn_idx.shape[0]), n_neighbors)
    col_indices = nn_idx.ravel().astype(int)
    weights = nn_dist.ravel()
    nn = sp_sparse.coo_matrix((weights, (row_indices, col_indices)), shape=(num_cells, num_cells))

    with subprocess.Popen(
        [CONVERT_BIN_PATH, "-i", "-", "-o", f"{LOUVAIN_PATH}/matrix"],
        stdin=subprocess.PIPE, text=True
    ) as proc:
        assert proc.stdin is not None
        for i, j in zip(nn.row, nn.col):
            proc.stdin.write(f"{i}\t{j}\n")
        proc.stdin.close()
        proc.wait()

    with open(f"{LOUVAIN_PATH}/out", 'w') as out_f:
        subprocess.call([
            LOUVAIN_BIN_PATH, f"{LOUVAIN_PATH}/matrix",
            "-q", "0", "-l", "-1", "-s", f"{random_seed}"
        ], stdout=out_f)

    labels = np.zeros(num_cells, dtype=np.int64)
    seen_idx = set()
    with open(f"{LOUVAIN_PATH}/out", 'r') as f:
        for line in f:
            used_bc_idx, cluster = map(int, line.strip().split(" "))
            if used_bc_idx in seen_idx:
                continue
            seen_idx.add(used_bc_idx)
            labels[used_bc_idx] = 1 + cluster

    pca_df = pd.DataFrame(pca_matrix)
    checked_cluster_pairs = set()

    while True:
        if len(np.unique(labels)) <= 1:
            break
        pca_df["cluster"] = labels
        medoids = pca_df.groupby("cluster").median().iloc[:, :-1].to_numpy()

        hc = linkage(medoids, method="complete")
        max_label = np.max(labels)
        any_merged = False

        for step in range(hc.shape[0]):
            left, right = hc[step, 0], hc[step, 1]
            if left <= max_label and right <= max_label:
                group0 = np.flatnonzero(labels == left)
                group1 = np.flatnonzero(labels == right)
                cluster_pair = tuple(sorted([frozenset(group0), frozenset(group1)]))
                if cluster_pair in checked_cluster_pairs:
                    continue
                checked_cluster_pairs.add(cluster_pair)

                sub_indices = np.concatenate((group0, group1))
                sub_matrix = sp_sparse.csc_matrix(expr_matrix[:, sub_indices])
                if sub_matrix.has_sorted_indices:
                    sub_matrix = sub_matrix.sorted_indices()
                params = compute_sseq_params_o3(sub_matrix.astype(np.int32), 0.995)
                group0_submat = np.arange(len(group0))
                group1_submat = np.arange(len(group0), len(group0) + len(group1))
                de_result = sseq_differential_expression_o3(
                    sub_matrix.astype(np.int32), group0_submat, group1_submat, params, 900
                )
                n_de_genes = np.sum(de_result["adjusted_p_values"] < 0.05)
                if n_de_genes == 0:
                    labels[labels == right] = left
                    labels[labels > right] -= 1
                    any_merged = True
                    break

        if not any_merged:
            break

    labels += 1
    relabel_order = np.argsort(np.argsort(-np.bincount(labels)))
    labels = 1 + relabel_order[labels]
    adata.obs[result_key] = pd.Categorical(labels.astype(str))


def preprocess_rna_data(
    data_list=None,
    marker=None,
    hk_genes_path="./He2020Nature_mouseHK.txt",
    min_genes=3,
    min_counts=3,
    exclude_prefixes=("Mt", "mt-", "Gm", "Rp"),
    exclude_rik=True
):
    if data_list is None:
        data_list = [
            "./RNA/B924-1.expr.bin_100.h5", "./RNA/B926-1.expr.bin_100.h5", "./RNA/B922-1.expr.bin_100.h5",
            "./RNA/B924-2.expr.bin_100.h5", "./RNA/B926-2.expr.bin_100.h5", "./RNA/B922-2.expr.bin_100.h5"
        ]
    if marker is None:
        marker = ["1-1", "2-1", "3-1", "1-2", "2-2", "3-2"]

    anndata_list = []
    for idx, data in enumerate(data_list):
        anndata_ = sc.read_h5ad(data)
        if marker is not None and idx < len(marker):
            anndata_.obs["marker"] = marker[idx]
        else:
            anndata_.obs["marker"] = f"marker_{idx}"
        anndata_list.append(anndata_)

    combined_data = ad.concat(anndata_list)
    combined_data.obs_names_make_unique()
    rnadata = combined_data.copy()

    for prefix in exclude_prefixes:
        rnadata.var[prefix] = rnadata.var_names.str.startswith(prefix)

    if exclude_rik:
        rnadata.var["rik"] = [True if "Rik" in str else False for str in rnadata.var_names]
    else:
        rnadata.var["rik"] = [False] * rnadata.shape[1]

    mask = np.ones(rnadata.shape[1], dtype=bool)
    for prefix in exclude_prefixes:
        mask &= ~rnadata.var[prefix].values
    mask &= ~rnadata.var["rik"].values
    fdata = rnadata[:, mask].copy()

    with open(hk_genes_path, 'r') as f:
        for line in f:
            hk_genes = line.strip().split('\t')
            break
    fdata = fdata[:, [gene for gene in fdata.var_names if gene not in hk_genes]]

    sc.pp.filter_cells(fdata, min_genes=min_genes)
    sc.pp.filter_cells(fdata, min_counts=min_counts)
    return fdata

def run_scvi_model(
    fdata: AnnData,
    model_save_path: str = "scvi.model",
    batch_key: str = "marker",
    n_hidden: int = 128,
    n_latent: int = 20,
    n_layers: int = 10,
    gene_likelihood: str = "poisson",
    latent_distribution: str = "normal",
    seed: int = 0
) -> scvi.model.SCVI:

    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(fdata, batch_key=batch_key)
    model = scvi.model.SCVI(
        fdata, 
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood=gene_likelihood,
        latent_distribution=latent_distribution
    )

    model.train(early_stopping=True, enable_progress_bar=True)

    model.save(model_save_path, save_anndata=True, overwrite=True)

    loaded_model = scvi.model.SCVI.load(model_save_path)
    return loaded_model

def run_cellcharter_clustering(
    fdata: AnnData,
    model,
    library_key: str = "marker",
    coord_type: str = "generic",
    spatial_key: str = "spatial",
    n_layers: int = 3,
    use_rep_scVI: str = "X_scVI",
    out_key: str = "X_cellcharter",
    sample_key: str = "time_point",
    n_clusters: tuple = (7, 20),
    max_runs: int = 10,
    random_state: int = 114514,
    cluster_obs_key: str = "cluster_cellcharter_autok"
) -> AnnData:

    fdata.obsm[use_rep_scVI] = model.get_latent_representation(fdata).astype(np.float32)
    sq.gr.spatial_neighbors(fdata, library_key=library_key, coord_type=coord_type, delaunay=True, spatial_key=spatial_key)
    cc.gr.remove_long_links(fdata)
    cc.gr.aggregate_neighbors(fdata, n_layers=n_layers, use_rep=use_rep_scVI, out_key=out_key, sample_key=sample_key)
    autok = cc.tl.ClusterAutoK(
        n_clusters=n_clusters,
        max_runs=max_runs,
        model_params=dict(random_state=random_state)
    )
    autok.fit(fdata, use_rep=out_key)
    fdata.obs[cluster_obs_key] = autok.predict(fdata, use_rep=out_key)
    return fdata

from typing import List, Optional, Tuple
from scipy.stats import pearsonr
from tqdm import tqdm

def nmf_consensus_programs(
    fdata: AnnData,
    n_top_genes: int = 10000,
    n_nmf_components: int = 100,
    nmf_max_iter: int = 2000,
    nmf_verbose: int = 1,
    min_clusters: int = 4,
    max_clusters: int = 20,
    n_resamples: int = 100,
    resample_frac: float = 0.8,
    program_defs: Optional[List[List[int]]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, object]:

    import consensusclustering as cc
    from sklearn.decomposition import NMF
    from sklearn.cluster import AgglomerativeClustering
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    
    sc.pp.highly_variable_genes(fdata, n_top_genes=n_top_genes, flavor='cell_ranger', inplace=True)
    hvg_data = fdata[:, fdata.var['highly_variable']].copy()
    
    cnt_matrix = hvg_data.X.toarray()
    model = NMF(n_components=n_nmf_components, verbose=nmf_verbose, max_iter=nmf_max_iter)
    cnt_matrix_trans = model.fit_transform(cnt_matrix)
    
    corr_matrix = np.zeros((cnt_matrix_trans.shape[1], cnt_matrix_trans.shape[1]))
    for i in tqdm(range(cnt_matrix_trans.shape[1])):
        for j in range(i, cnt_matrix_trans.shape[1]):
            corr_matrix[i, j] = pearsonr(cnt_matrix_trans[:, i], cnt_matrix_trans[:, j])[0]
    for i in tqdm(range(corr_matrix.shape[1])):
        for j in range(i):
            corr_matrix[i, j] = corr_matrix[j, i]
    
    cc_model = cc.ConsensusClustering(
        AgglomerativeClustering(),
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        n_resamples=n_resamples,
        resample_frac=resample_frac,
        k_param='n_clusters'
    )
    cc_model.fit(corr_matrix)
    best_k = cc_model.best_k()

    if program_defs is None:
        dist_matrix = 1 - corr_matrix
        np.fill_diagonal(dist_matrix, 0)
        condensed_dist = squareform(dist_matrix, checks=False)
        Z = linkage(condensed_dist, method='average')
        clusters = fcluster(Z, t=best_k, criterion='maxclust')
        program_defs = []
        for k in range(1, best_k + 1):
            indices = np.where(clusters == k)[0].tolist()
            if indices:
                program_defs.append(indices)

    combined_matrix_trans = np.zeros((len(program_defs), len(fdata.obs_names)))
    for i, program in enumerate(program_defs):
        combined_matrix_trans[i, :] = cnt_matrix_trans[:, program].sum(axis=1)
        combined_matrix_trans[i, :] = combined_matrix_trans[i, :] / combined_matrix_trans[i, :].sum() * 1e3

    return combined_matrix_trans, cnt_matrix_trans, corr_matrix, cc_model

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
sc.settings.verbosity = 3

def calculate_composition(
    adata: AnnData,
    groupby_key: str,
    category_key: str
) -> pd.DataFrame:
    """
    Calculate the composition of categories within each group.

    Parameters:
    adata: AnnData object
    groupby_key: Key in adata.obs to group by (e.g., 'louvain')
    category_key: Key in adata.obs for categories (e.g., 'annotations')

    Returns:
    Pandas DataFrame with composition data.
    """
    # Create a DataFrame from the AnnData object
    data = pd.DataFrame(adata.obs)
    # Calculate composition
    composition = data.groupby(groupby_key)[category_key].value_counts(normalize=True).unstack(fill_value=0)
    return composition


def plot_composition(
    composition_df: pd.DataFrame,
    title: str,
    x: str,
    y: str,
    palette: str = 'viridis'
) -> None:
    """
    Plot the composition data as a stacked bar plot.
    Parameters:
    - composition_df: DataFrame with composition data
    - title: Title for the plot
    - x: X-axis label
    - y: Y-axis label
    - palette: Color palette
    Return:
    - None
    """
    if isinstance(palette, str):
        cmap = plt.get_cmap(palette)
        colors = cmap(np.linspace(0, 1, composition_df.shape[1]))
    elif isinstance(palette, list):
        colors = palette
    else:
        raise ValueError("Palette should be a string (colormap name) or a list of colors")
    
    ax = composition_df.plot(kind='bar', stacked=True, figsize=(6, 4), color=colors, fontsize=8)
    
    plt.title(title)
    plt.ylabel('Proportion')
    plt.xlabel(x)

    plt.legend(title=y, bbox_to_anchor=(1.05, 1), loc='upper left')

    sns.despine()

    ax.yaxis.grid(False)
    ax.xaxis.grid(False)

    plt.tight_layout()
    plt.show()

def compute_relative_abundance(
    adata: AnnData,
    dict_use: dict[str, list[str]]
) -> dict[str, float]:
    relative_abundances: dict[str, float] = {}
    for cell_type, markers in tqdm(dict_use.items()):
        avg_expressions = []
        for marker in markers:
            if marker not in adata.var_names:
                continue
            expressed_values = adata[:, marker].X[adata[:, marker].X > 0]
            if isinstance(expressed_values, np.ndarray):
                num_nonzero = len(expressed_values)
            else:
                num_nonzero = expressed_values.getnnz()
            if num_nonzero > 0:
                avg_expressions.append(np.mean(expressed_values))
            else:
                avg_expressions.append(0)
        relative_abundance = np.log10(sum(avg_expressions) + 1e-10)
        relative_abundances[cell_type] = 1 / relative_abundance
    return relative_abundances


def annotate_cells_stage(
    marker_dict_use: dict[str, list[str]],
    adata: AnnData,
    cells_to_annotate: Optional[list[str]] = None
) -> tuple[dict[str, str], list[str]]:
    WCT = compute_relative_abundance(adata, marker_dict_use)
    annotations: dict[str, str] = {}
    unannotated_cells: list[str] = []

    if cells_to_annotate is None:
        cells_to_annotate: list[str] = adata.obs_names
        data_submatrix = adata.X
    else:
        data_submatrix = adata[cells_to_annotate, :].X

    all_scores = np.zeros((len(cells_to_annotate), len(marker_dict_use)))
    for idx, (cell_type, markers) in enumerate(marker_dict_use.items()):
        valid_markers_indices = [
            adata.var_names.get_loc(marker)
            for marker in markers
            if marker in adata.var_names
        ]
        marker_matrix = data_submatrix[:, valid_markers_indices]
        presence_matrix = (marker_matrix > 0).astype(int)
        scores = presence_matrix.sum(axis=1) * WCT[cell_type]

        all_scores[:, idx] = scores.ravel()
    max_scores = np.max(all_scores, axis=1)
    max_score_indices = np.argmax(all_scores, axis=1)
    cell_types = list(marker_dict_use.keys())
    annotations = np.array(
        [
            cell_types[idx] if score > 0 else "Others"
            for idx, score in zip(max_score_indices, max_scores)
        ]
    )
    unannotated = [cells_to_annotate[i] for i in np.where(annotations == "Others")[0]]
    return dict(zip(cells_to_annotate, annotations)), unannotated

def export_tpm_by_clone(
    adata_path: str = "/data200T/SPACseq/HD/output/clones_with_labels.h5ad",
    output_path: str = "./tpm_by_clone.txt"
):
    adata = sc.read_h5ad(adata_path)
    counts = adata.X
    if not isinstance(counts, np.ndarray):
        counts = counts.toarray()

    adata.var_names = adata.var_names.str.upper()
    adata.obs["clone"] = adata.obs["clone"].astype(str)
    clones = adata.obs["clone"].unique()

    tpm_df = pd.DataFrame(index=adata.var_names, columns=clones)
    for clone in clones:
        clone_indices = adata.obs.index[adata.obs["clone"] == clone]
        clone_counts = counts[adata.obs.index.isin(clone_indices), :].sum(axis=0)
        tpm = (clone_counts / clone_counts.sum()) * 1e6
        tpm_df[clone] = tpm

    tpm_df.index.name = "Gene"
    tpm_df.index.name = None
    tpm_df.to_csv(output_path, sep="\t")
