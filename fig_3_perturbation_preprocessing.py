import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from anndata import AnnData

from A_utils import remove_mito_ribo_hk_lnc_genes, combine_guide_replicates

def prepare_guide_data(
    samples: tuple[str, ...] = ('B922-1', 'B922-2', 'B924-1', 'B924-2', 'B926-1', 'B926-2'),
    h5_dir: str = './h5_files/',
):
    adata = {}
    for sample in samples:
        adata[sample] = sc.read_h5ad(f'{h5_dir}{sample}.bin200.h5')
        adata[sample].obs['marker'] = sample
    adata_combined = ad.concat(adata.values(), join='outer')
    adata_combined.obs_names_make_unique()

    fdata = adata_combined[:, ~adata_combined.var_names.str.startswith('sg')].copy()
    gdata = adata_combined[:, adata_combined.var_names.str.startswith('sg')].copy()

    # Merge sgNTC and sgnon-targeting genes, as they are named differently in the data
    gdata.X[:, gdata.var_names == 'sgNTC'] = (
        gdata[:, 'sgNTC'].X.toarray() +
        gdata[:, 'sgnon-targeting_1_gene'].X.toarray() +
        gdata[:, 'sgnon-targeting_2_gene'].X.toarray()
    )
    gdata = gdata[:, ~gdata.var_names.str.startswith('sgnon-targeting')].copy()

    gdata = combine_guide_replicates(gdata)

    fdata = remove_mito_ribo_hk_lnc_genes(fdata)
    fdata.layers['counts'] = fdata.X.copy()

    sc.pp.normalize_total(fdata, inplace=True, target_sum=1e4)
    sc.pp.log1p(fdata)

    return fdata, gdata

# Align RNA and guide data
def clean_guide_data(
    rnadata: AnnData,
    guidedata: AnnData
) -> tuple[AnnData, AnnData]:

    rnadata.obs["cov"] = [str(marker) + '_' + str(array[0]) + "-" + str(array[1]) for marker, array in zip(rnadata.obs["marker"], rnadata.obsm["spatial"])]
    guidedata.obs["cov"] = [str(marker) + '_' + str(array[0]) + "-" + str(array[1]) for marker, array in zip(guidedata.obs["marker"], guidedata.obsm["spatial"])]

    common_cov = np.intersect1d(rnadata.obs['cov'], guidedata.obs['cov'])
    guidedata.obs_names = guidedata.obs['cov']
    rnadata.obs_names = rnadata.obs['cov']

    guidedata_filtered = guidedata[common_cov].copy()
    sorted_obs_names = rnadata.obs_names.sort_values()
    guidedata_filtered = guidedata_filtered[sorted_obs_names].copy()
    rnadata = rnadata[sorted_obs_names].copy()
    missing_cov = rnadata.obs['cov'][~rnadata.obs['cov'].isin(common_cov)]
    if len(missing_cov) > 0:
        missing_obs = pd.DataFrame({'marker': [cov.split('_')[0] for cov in missing_cov]}, index=missing_cov)
        missing_obsm_spatial = np.array([[float(cov.split('_')[1].split('-')[0]), float(cov.split('_')[1].split('-')[1])] for cov in missing_cov])
        missing_X = np.zeros((len(missing_cov), guidedata.shape[1]))

        missing_guidedata = ad.AnnData(X=missing_X, obs=missing_obs, obsm={'spatial': missing_obsm_spatial})
        guidedata_filtered = guidedata_filtered.concatenate(missing_guidedata)
    return rnadata, guidedata_filtered

