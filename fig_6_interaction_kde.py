import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm
from scipy.ndimage import gaussian_filter1d
from scipy.stats import ttest_ind

def kde_interaction_analysis(
    rna_h5ad='./RNA/reclustered.h5',
    spatial_h5ad='./RNA/spatial.combined.bin100.h5',
    marker='A',
    y_threshold=4000,
    gene_expr='Cxcl12',
    guide1='sgBhlhe40',
    guide2='sgNon-targeting',
    grid_size=50,
    smooth_sigma=18
):
    """
    Performs KDE-based interaction analysis on spatial transcriptomics data.

    Parameters:
        rna_h5ad (str): Path to single-cell RNA h5ad file.
        spatial_h5ad (str): Path to spatial data h5ad file.
        marker (str): Value in 'marker' column to select cells/spots.
        y_threshold (int): Minimum value for spatial y to subset data.
        gene_expr (str): Gene name for KDE expression.
        guide1, guide2 (str): Guide column names in spatial data to compare.
        grid_size (int): Number of grid points for KDE mesh.
        smooth_sigma (int): Smoothing parameter for gaussian_filter1d.

    Returns:
        df (pd.DataFrame): DataFrame of smoothed and normalized guide abundances per KDE level.
        t_stat (float): t-statistic from t-test.
        p_value (float): p-value from t-test.
    """

    fdata = sc.read_h5ad(rna_h5ad)
    gdata = sc.read_h5ad(spatial_h5ad)
    gdata = gdata[:, gdata.var_names.str.startswith('sg')].copy()

    # Sample 1
    plot_data = fdata[(fdata.obs['marker'] == marker) & (fdata.obsm['spatial'][:, 1] > y_threshold)]
    plot_guide_data = gdata[(gdata.obs['marker'] == marker) & (gdata.obsm['spatial'][:, 1] > y_threshold)]

    # KDE on X/Y and expression
    kde = sm.nonparametric.KDEMultivariate(
        data=[
            plot_data.obsm['spatial'][:, 0],
            plot_data.obsm['spatial'][:, 1],
            plot_data[:, gene_expr].X.toarray().flatten()
        ],
        var_type='ccc',
        bw='normal_reference'
    )

    # Generate grid points
    x_grid, y_grid, c_grid = np.meshgrid(
        np.linspace(plot_data.obsm['spatial'][:, 0].min(), plot_data.obsm['spatial'][:, 0].max(), grid_size),
        np.linspace(plot_data.obsm['spatial'][:, 1].min(), plot_data.obsm['spatial'][:, 1].max(), grid_size),
        np.linspace(
            plot_data[:, gene_expr].X.toarray().flatten().min(),
            plot_data[:, gene_expr].X.toarray().flatten().max(),
            grid_size
        )
    )
    grid_coords = np.column_stack([x_grid.ravel(), y_grid.ravel(), c_grid.ravel()])

    # KDE values and sum across expression
    kde_values = kde.pdf(grid_coords).reshape(x_grid.shape)
    kde_values = kde_values.sum(axis=2)

    threshold = np.quantile(kde_values, 0.95)

    kde_df = pd.DataFrame({
        'x': x_grid[:, :, 0].flatten() // grid_size,
        'y': y_grid[:, :, 0].flatten() // grid_size,
        'kde': kde_values.flatten()
    })
    kde_df = kde_df.groupby(['x', 'y']).agg({'kde': 'mean'}).reset_index()
    kde_df['level'] = kde_df['kde'].rank() // 10
    kde_df['level'] = kde_df['level'].astype(int)

    guide_df = pd.DataFrame({
        'x': plot_guide_data.obsm['spatial'][:, 0] // grid_size,
        'y': plot_guide_data.obsm['spatial'][:, 1] // grid_size,
        guide1: plot_guide_data[:, guide1].X.toarray().flatten(),
        guide2: plot_guide_data[:, guide2].X.toarray().flatten()
    })
    guide_df = guide_df.groupby(['x', 'y']).agg({guide1: 'sum', guide2: 'sum'}).reset_index()
    df = pd.merge(kde_df, guide_df, on=['x', 'y'], how='left')
    df = df.dropna().groupby('level').agg({'kde': 'mean', guide1: 'sum', guide2: 'sum'}).reset_index()
    df[guide1] = df[guide1] / df[guide1].sum()
    df[guide2] = df[guide2] / df[guide2].sum()

    # Smooth data
    df[guide1] = gaussian_filter1d(df[guide1], sigma=smooth_sigma)
    df[guide2] = gaussian_filter1d(df[guide2], sigma=smooth_sigma)

    t_stat, p_value = ttest_ind(df[df.kde > threshold][guide1], df[df.kde > threshold][guide2], equal_var=False)

    return df, t_stat, p_value