import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

def compute_log_normalized_guide_cluster_matrix(
    guide_data: ad.AnnData,
    non_targeting_label: str = 'sgnon-targeting'
) -> pd.DataFrame:
    """
    Computes a log-normalized guide-cluster matrix.

    Parameters:
    - guide_data: AnnData object with guide counts and 'cluster' in .obs
    - non_targeting_label: name of the non-targeting guide (row)
    Return:
    - p_df: pandas DataFrame, log10-normalized and ratioed to non-targeting
    """
    c_df = pd.DataFrame(guide_data.X, columns=guide_data.var_names)
    c_df["cluster"] = guide_data.obs["cluster"].tolist()
    c_df = c_df.groupby(["cluster"]).sum()

    p_df = c_df.T.div(c_df.T.sum(axis=1), axis=0)
    p_df = p_df.div(p_df.loc[non_targeting_label, :], axis=1)
    p_df = np.log10(p_df)

    return p_df

from scipy.stats import chi2_contingency

def compute_cluster_guide_stats(
    mdf: pd.DataFrame,
    gdata: ad.AnnData,
    marker: str = '2-1',
    non_targeting_label: str = 'sgnon-targeting'
):
    """
    Compute Chi-square p-values and log2 fold-change (l2fc) for each cluster-guide pair,
    and annotate p-values for significance.

    Parameters:
        mdf: DataFrame with columns ['marker', 'cluster', 'variable', 'value']
        gdata: object providing .var_names for guide identifiers
        marker: string marker to use for group inclusion
        non_targeting_label: string guide name for non-targeting control

    Returns:
        chi_df: DataFrame of p-values (clusters x guides)
        pdf: DataFrame of log2 fold-changes (clusters x guides)
        annot_df: DataFrame of string annotations mapping p-values to asterisks
    """
    p_dict = {}

    for cluster in mdf.cluster.unique():
        p_dict[cluster] = []
        for guide in mdf.variable.unique():
            if guide == non_targeting_label:
                continue
            int_df = mdf[
                (mdf.marker == marker) &
                (mdf.cluster == cluster) &
                (np.isin(mdf.variable, [guide, non_targeting_label]))
            ][['variable', 'value']].set_index('variable')
            res_df = mdf[
                (mdf.marker != marker) &
                (mdf.cluster != cluster) &
                (np.isin(mdf.variable, [guide, non_targeting_label]))
            ].set_index('cluster')[['variable', 'value']].groupby('variable').sum()
            chi_df_inner = pd.concat([int_df, res_df], axis=1)
            chi_df_inner.columns = [f'cluster_{cluster}', 'rest']
            # If there are missing guides (rare), pad with zeros
            for idx in [guide, non_targeting_label]:
                if idx not in chi_df_inner.index:
                    chi_df_inner.loc[idx] = [0, 0]
            chi_df_inner = chi_df_inner.loc[[guide, non_targeting_label]]  # Ensure row order
            try:
                pval = chi2_contingency(chi_df_inner)[1]
            except:
                pval = np.nan
            p_dict[cluster].append(pval)

    guides = [guide for guide in gdata.var_names if guide != non_targeting_label]
    chi_df = pd.DataFrame(p_dict, index=guides).T

    pdf = mdf.set_index(['cluster', 'variable'])['value'].unstack()
    pdf = pdf.div(pdf.sum(axis=0), axis=1)
    pdf = pdf.div(pdf[non_targeting_label], axis=0).loc[:, pdf.columns != non_targeting_label]

    plot_df = pd.DataFrame({
        'p-value':  chi_df.melt().value.tolist(),
        'FC': list(np.log2(pdf.melt().value + 1e-9))
    }, index=[str(a) + ' c' + str(b) for a, b in zip(chi_df.melt().variable, chi_df.melt(ignore_index=False).index)])

    plot_df = plot_df[plot_df.FC > np.log2(1e-6)]

    plot_df.sort_values(by='p-value', ascending=True, inplace=True)
    plot_df['rank'] = np.arange(len(plot_df))
    plot_df['padj'] = plot_df['p-value'] * len(plot_df) / (plot_df['rank'] + 1)
    plot_df['padj'] = -np.log10(plot_df['padj'].clip(upper=1))

    return plot_df

def compute_degs_and_l2fc_pvals(
    annotation,
    annotation_guide,
    cdata,
    output_deg_csv_path: str = './deg.csv',
    pval_cutoff: float = 1,
    top_n: int = 300000
):
    """
    Process single-cell data to compute DEGs, log2 fold changes, and p-value DataFrames.

    Args:
        annotation: AnnData object containing cell annotations.
        annotation_guide: AnnData or view-like structure with columns for sgRNAs.
        cdata: AnnData object for Scanpy DEG analysis.
        output_deg_csv_path: Output file path for DEG results.
        pval_cutoff: p-value threshold for filtering DEGs.
        top_n: Number of top genes to retain per group.

    Returns:
        p_df: DataFrame containing adjusted p-values for Control and sgCd44 groups.
    """
    # Subset groups by guide identity
    guide_data = annotation[(annotation_guide[:, 'sgCd44'].X > 0) & (annotation_guide[:, 'sgnon-targeting'].X == 0) & (annotation_guide[:, 'sgGata3'].X == 0)].copy()
    ntc_data = annotation[(annotation_guide[:, 'sgCd44'].X == 0) & (annotation_guide[:, 'sgnon-targeting'].X > 0) & (annotation_guide[:, 'sgGata3'].X == 0)].copy()
    pos_data = annotation[(annotation_guide[:, 'sgCd44'].X == 0) & (annotation_guide[:, 'sgnon-targeting'].X == 0) & (annotation_guide[:, 'sgGata3'].X > 0)].copy()

    guide_data.obs['guide'] = 'sgCd44'
    ntc_data.obs['guide'] = 'Control'
    pos_data.obs['guide'] = 'sgGata3'

    sc.tl.rank_genes_groups(guide_data, groupby='guide', method='t-test')
    sc.tl.rank_genes_groups(ntc_data, groupby='guide', method='t-test')
    sc.tl.rank_genes_groups(pos_data, groupby='guide', method='t-test')

    result = guide_data.uns['rank_genes_groups']
    result = ntc_data.uns['rank_genes_groups']
    result = pos_data.uns['rank_genes_groups']

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

        filtered_data = data[data['padj'] < pval_cutoff]
        sorted_data = filtered_data.sort_values(by='score', ascending=False)
        top_genes = sorted_data.head(top_n)
        top_genes_df = pd.concat([top_genes_df, top_genes], ignore_index=True)

    deg = top_genes_df["Gene"].unique().tolist()
    top_genes_df.to_csv(output_deg_csv_path, index=None, sep='\t')
    # Optionally inspect the Spp1 result, but this has no effect if not used
    _ = top_genes_df[top_genes_df.Gene == 'Spp1']

    # Select genes for control and guide for visual comparison
    genes = np.unique(
        sc.get.rank_genes_groups_df(cdata, group='sgCd44').head(12).names.tolist()
        + sc.get.rank_genes_groups_df(cdata, group='Control').head(12).names.tolist()
    )

    # Log2 fold change DataFrame
    ntc_df = sc.get.rank_genes_groups_df(cdata, group='Control') \
        .loc[np.isin(sc.get.rank_genes_groups_df(cdata, group='Control').names, genes), :][['names', 'logfoldchanges']] \
        .set_index('names')
    ntc_df.columns = ['Control']
    guide_df = sc.get.rank_genes_groups_df(cdata, group='sgCd44') \
        .loc[np.isin(sc.get.rank_genes_groups_df(cdata, group='sgCd44').names, genes), :][['names', 'logfoldchanges']] \
        .set_index('names')
    guide_df.columns = ['sgCd44']
    l2fc_df = pd.concat([ntc_df, guide_df], axis=1)

    # Adjusted p-value DataFrame
    ntc_df = sc.get.rank_genes_groups_df(cdata, group='Control') \
        .loc[np.isin(sc.get.rank_genes_groups_df(cdata, group='Control').names, genes), :][['names', 'pvals_adj']] \
        .set_index('names')
    ntc_df.columns = ['Control']
    guide_df = sc.get.rank_genes_groups_df(cdata, group='sgCd44') \
        .loc[np.isin(sc.get.rank_genes_groups_df(cdata, group='sgCd44').names, genes), :][['names', 'pvals_adj']] \
        .set_index('names')
    guide_df.columns = ['sgCd44']
    p_df = pd.concat([ntc_df, guide_df], axis=1)

    return p_df

def compute_unified_levels(
    x_coords,
    y_coords,
    weights_guide,
    weights_control,
    bw_adjust: float = 0.4,
    num_levels: int = 7,
    level_min: float = 1e-9,
    level_max: float = 1e-8
) -> tuple[np.ndarray, np.ndarray]:
    linear_mapping = np.linspace(level_min, level_max, num_levels)

    sgCd44_quantiles = []
    sgNTC_quantiles = []

    def get_density_range(x, y, weights, bw_adjust):
        from scipy.stats import gaussian_kde
        # Ensure float type for arithmetic to avoid dtype casting errors
        xmin, xmax = np.min(x).astype(float), np.max(x).astype(float)
        ymin, ymax = np.min(y).astype(float), np.max(y).astype(float)
        x_range = xmax - xmin
        y_range = ymax - ymin
        xmin = xmin - x_range * 0.05
        xmax = xmax + x_range * 0.05
        ymin = ymin - y_range * 0.05
        ymax = ymax + y_range * 0.05

        xx, yy = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]
        positions = np.vstack([xx.ravel(), yy.ravel()])

        kde = gaussian_kde(np.vstack([x, y]), weights=weights)
        kde.set_bandwidth(kde.factor * bw_adjust)
        density = kde(positions)

        return density

    densities_sgGuide = get_density_range(x_coords, y_coords, weights_guide, bw_adjust)
    densities_sgControl = get_density_range(x_coords, y_coords, weights_control, bw_adjust)

    quantiles = np.linspace(0, 1, 1001)[1:]  # 1% to 100%
    sgGuide_values = np.quantile(densities_sgGuide.flatten(), quantiles)
    sgControl_values = np.quantile(densities_sgControl.flatten(), quantiles)

    sgGuide_quantiles = []
    sgControl_quantiles = []
    for target in linear_mapping:
        nearest_idx = np.abs(sgGuide_values - target).argmin()
        nearest_quantile = quantiles[nearest_idx]
        nearest_value = sgGuide_values[nearest_idx]
        sgCd44_quantiles.append(nearest_quantile)

        nearest_idx = np.abs(sgControl_values - target).argmin()
        nearest_quantile = quantiles[nearest_idx]
        nearest_value = sgControl_values[nearest_idx]
        sgNTC_quantiles.append(nearest_quantile)

    sgCd44_quantiles = np.array(sgCd44_quantiles + [1.0])
    sgNTC_quantiles = np.array(sgNTC_quantiles + [1.0])

    return sgCd44_quantiles, sgNTC_quantiles

def process_samples_for_quantiles(
    fdata, fdata_bin200, gdata_raw, region, guide_name,
    samples: list[str] = None, bw_adjust: float = 0.4, num_levels: int = 7, 
    level_min: float = 1e-9, level_max: float = 1e-8
):
    """
    For each sample, extract spatial coordinates and guide weights,
    and compute unified quantile levels.
    
    Args:
        fdata: AnnData object with cell metadata and spatial info
        fdata_bin200: AnnData object (not used in this function)
        gdata_raw: AnnData object with guide data
        region: ((xmin, xmax), (ymin, ymax)), spatial region
        samples: list of sample names (if None, uses default set)
        bw_adjust, num_levels, level_min, level_max: params for compute_unified_levels
        
    Returns:
        results: dict mapping sample name to (sgCd44_quantiles, sgNTC_quantiles)
    """
    if samples is None:
        samples = ['1-1', '1-2', '2-1', '2-2', '3-1', '3-2']

    results = {}
    for SAMPLE in samples:
        sample_idx = np.where(fdata_bin200.obs['marker'].values == f'{SAMPLE}')[0]

        # 提取坐标和weights
        plot_data = fdata[(fdata.obs['marker'] == SAMPLE) & 
                          (fdata.obsm['spatial'][:, 0] >= region[0][0]) &
                          (fdata.obsm['spatial'][:, 0] <= region[0][1]) &
                          (fdata.obsm['spatial'][:, 1] >= region[1][0]) &
                          (fdata.obsm['spatial'][:, 1] <= region[1][1])]
        plot_guide_data = gdata_raw[(gdata_raw.obs['marker'] == SAMPLE) & 
                                    (gdata_raw.obsm['spatial'][:, 0] >= region[0][0]) &
                                    (gdata_raw.obsm['spatial'][:, 0] <= region[0][1]) &
                                    (gdata_raw.obsm['spatial'][:, 1] >= region[1][0]) &
                                    (gdata_raw.obsm['spatial'][:, 1] <= region[1][1])]
        x_coords = plot_data.obsm['spatial'][:, 0].astype(float)  # Ensure float dtype
        y_coords = plot_data.obsm['spatial'][:, 1].astype(float)
        weights_cd44 = plot_guide_data[:, guide_name].X.toarray().flatten()
        weights_ntc = plot_guide_data[:, 'sgnon-targeting'].X.toarray().flatten()
        sgCd44_quantiles, sgNTC_quantiles = compute_unified_levels(
            x_coords, y_coords, weights_cd44, weights_ntc, 
            bw_adjust, num_levels, level_min, level_max
        )
        results[SAMPLE] = (sgCd44_quantiles, sgNTC_quantiles)
    return results

def compute_cluster_timepoint_proportions(
    fdata: ad.AnnData,
    cluster_col: str = 'cluster_cellcharter_givenk',
    timepoint_col: str = 'timepoint'
) -> pd.DataFrame:
    """
    Computes the proportion of bins for each cluster at each timepoint.

    Args:
        fdata: AnnData object with .obs containing cluster and timepoint columns
        cluster_col: Name of column in .obs denoting cluster assignments
        timepoint_col: Name of column in .obs denoting timepoint

    Returns:
        prop_df: DataFrame with percent composition of each cluster at each timepoint,
            and a column for re-labeled spatial niche labels.
    """
    cnt_df = fdata.obs.groupby([cluster_col, timepoint_col]).size().unstack()
    prop_df = cnt_df.div(cnt_df.sum(axis=0), axis=1) * 100  # percentage
    prop_df = prop_df.reset_index().melt(id_vars=cluster_col, var_name='Time', value_name='Proportion (%)')

    prop_df = prop_df.sort_values(by=[cluster_col, 'Time'])
    prop_df['Spatial niches'] = prop_df[cluster_col].astype(int).map(lambda x: f"N{x+1}")
    return prop_df