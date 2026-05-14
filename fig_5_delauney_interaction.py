import numpy as np
import pandas as pd
import squidpy as sq
from anndata import AnnData

def compute_binary_enrichment(
    inter_data: AnnData,
    annotation_: AnnData,
    field_name: str = "",
    radius: int = 25
) -> pd.DataFrame:
    """
    Computes the binary enrichment matrix using neighborhood enrichment for a given field.
    
    Parameters:
        - inter_data: AnnData object containing data for interaction analysis.
        - annotation_: AnnData object containing annotation, used for spatial neighbors.
        - field_name: The field/column in .obs to use for grouping.
        - radius: Radius for spatial neighbors computation.
    
    Returns:
        - binary_enrichment_df: DataFrame representing the binary enrichment matrix.
    """
    inter_data.obs[field_name] = inter_data.obs[field_name].astype('category')
    if "spatial_neighbors" not in annotation_.uns.keys():
        sq.gr.spatial_neighbors(annotation_, delaunay=False, radius=radius, coord_type="generic")
    sq.gr.interaction_matrix(inter_data, field_name)
    sq.gr.nhood_enrichment(inter_data, field_name)

    enrichment_matrix = inter_data.uns[f"{field_name}_nhood_enrichment"]['count']
    zscores = inter_data.uns[f"{field_name}_nhood_enrichment"]['zscore']
    relative_enrichment = enrichment_matrix > 0
    significant_enrichment = (zscores > np.quantile(zscores, 0.5)) & relative_enrichment

    binary_enrichment_matrix = np.where(significant_enrichment, 1, 0)
    cluster_names = inter_data.obs[f"{field_name}"].cat.categories
    binary_enrichment_df = pd.DataFrame(binary_enrichment_matrix, 
                                        index=cluster_names, 
                                        columns=cluster_names)
    return binary_enrichment_df