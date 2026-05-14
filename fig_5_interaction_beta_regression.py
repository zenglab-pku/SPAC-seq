import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV
import seaborn as sns

markers = {
    "Neutrophils": ["PTPRC", "S100A9", "MPO", "CXCR2", "CSF3R", "AQP9"],
    "Macrophages": ["PTPRC", "CD14", "ITGAM", "ITGAX", "TREM2", "FCGR3A", "FCGR3B", "CCR2", "CD163", "CD68", "ARG1", "ADGRE1", "CD33", "CD68", "C1QA", "CD163",  "C1QB", "C1QC",  "MARCO", "CD81"],
    "DC": ["PTPRC", "ITGAM", "ITGAX", "CD14", "CD33", "CD68", "CD1C", "CD209", "CCL22", "CD86", "CLEC9A", "CLEC10A", "XCR1", "LAMP3"],
    "Endothelial": ["VWF", "PECAM1", "CLDN5","CDH5","MGP"],
    "Epithelial": ["CDH1", "EPCAM", "CDKN2A", "KRT5", "KRT17", "MUC5B"],
    "Neutrophils": ["PTPRC", "S100A9", "MPO", "CXCR2", "CSF3R", "AQP9"],
    "Malignant": ["MUC1", "FXYD3", "MKI67", "MDK", "SOX4", "S100P",  "EPCAM", "ATP1B1"],
    "Monocyte": ["FCN1", "LYZ", "CD14", "S100A8", "CFP"],
    "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1","MMP2", "MYL9", "DCN","FAP"],
}

def compute_marker_beta_regression(
    fdata_path: str,
    cdata_path: str,
    markers: dict,
    cluster_value: int = 3,
    gene_mouse_fn=None,
    alpha_list=[0.001, 0.005, 0.01, 0.05, 0.1],
    max_iter: int = 10000,
    l1_ratio: float = 0.5,
    verbose: bool = True
):
    """
    Compute marker set scores and perform beta regression between RNA and guide data for a specific cluster.

    Args:
        fdata_path (str): Path to clustered AnnData file (RNA bins).
        cdata_path (str): Path to guide AnnData file (guide bins).
        markers (dict): Dictionary mapping marker names to gene symbol lists.
        cluster_value (int, optional): Cluster value (from 'cluster_cellcharter_givenk') to filter by. Default is 3.
        gene_mouse_fn (callable, optional): Function to map marker gene list to mouse gene symbols.
        alpha_list (list, optional): List of alpha values for ElasticNet regularization.
        max_iter (int, optional): Max iterations for ElasticNet.
        l1_ratio (float, optional): L1 ratio for ElasticNet.
        verbose (bool, optional): Verbosity flag.

    Returns:
        pd.DataFrame: Regression coefficients dataframe (guides × marker_score terms).
    """
    fdata = ad.read_h5ad(fdata_path)
    cdata = ad.read_h5ad(cdata_path)

    for marker in markers:
        gene_list = markers[marker]
        genes = gene_mouse_fn(gene_list) if gene_mouse_fn else gene_list
        mean_expression = np.zeros(len(fdata.obs_names))
        for gene in genes:
            if gene in fdata.var_names:
                mean_expression += fdata[:, gene].X.toarray().flatten()
        fdata.obs[f"{marker}_score"] = mean_expression

    mask = fdata.obs["cluster_cellcharter_givenk"] == cluster_value

    def beta_regression(
        rnadata, guidedata, expr_term=[], expr_layer=None, 
        alpha_list=[0.001, 0.005, 0.01, 0.05, 0.1], max_iter=10000, l1_ratio=0.5, verbose=True
    ):
        """
        Perform ElasticNet regression for guide vs. bin marker scores or gene expression.

        Args:
            rnadata (AnnData): AnnData object containing RNA info.
            guidedata (AnnData): AnnData object with guide abundance.
            expr_term (list): Terms in .obs or .var_names to use as features (objectives).
            expr_layer (str, optional): Data matrix layer to use.
            alpha_list (list): List of alpha parameter values for ElasticNet.
            max_iter (int): Max iterations for ElasticNet.
            l1_ratio (float): L1/L2 penalty ratio.
            verbose (bool): Verbosity flag.
        Returns:
            pd.DataFrame: Regression coefficient matrix [guides × objectives].
        """
        if len(rnadata.obs_names) != len(guidedata.obs_names):
            print("Error!, data bin count doesn't match!")
            return None
        if type(guidedata.X) == "scipy.sparse._csr.csr_matrix":
            guide_matrix = guidedata.X.toarray()
        else:
            guide_matrix = guidedata.X

        # Inputs are in .obs (marker scores)
        if all([obj in rnadata.obs for obj in expr_term]):
            if verbose: print("Extracting bin metric")
            score_df = rnadata.obs[expr_term]
            try:
                score_matrix = np.array(score_df.astype(int))
            except:
                print("Error converting obs to np.float matrix!")
                return None
        # Inputs are genes
        elif all([obj in rnadata.var_names for obj in expr_term]):
            if verbose: print("Extracting expression metric")
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
        if verbose: 
            print(
                f"Regressing with {guide_matrix.shape[0]} samples, "
                f"{guide_matrix.shape[1]} guides, {score_matrix.shape[1]} objectives."
            )
        param_grid = { 'alpha': alpha_list }
        for i in tqdm(range(guide_matrix.shape[1])):
            guide_vector = guide_matrix[:, i]
            regr = ElasticNet(max_iter=max_iter, l1_ratio=l1_ratio)
            grid_regr = GridSearchCV(regr, param_grid, scoring='neg_mean_squared_error', cv=5, n_jobs=-1, verbose=0)
            grid_regr.fit(X=score_matrix, y=guide_vector)
            regr = grid_regr.best_estimator_
            regr.fit(X=score_matrix, y=guide_vector)
            corr_matrix[i] = regr.coef_
        corr_df = pd.DataFrame(corr_matrix, index=guidedata.var_names, columns=expr_term)
        if verbose: print("Regression done!")
        return corr_df

    marker_score_terms = [f"{marker}_score" for marker in markers]
    r_df = beta_regression(fdata[mask], cdata[mask], expr_term=marker_score_terms, 
                           alpha_list=alpha_list, max_iter=max_iter, l1_ratio=l1_ratio, verbose=verbose)
    return r_df