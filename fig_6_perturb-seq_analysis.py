import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import pandas as pd

from scipy.cluster.hierarchy import fcluster

from A_utils import remove_mito_ribo_hk_lnc_genes
def perturb_data_preprocessing():
    adata = sc.read_10x_mtx(
        './B646/filtered_feature_bc_matrix/',
        var_names='gene_symbols',
        cache=True,
        gex_only=False
    )
    bdata = sc.read_10x_mtx(
        './B647/filtered_feature_bc_matrix/',
        var_names='gene_symbols',
        cache=True,
        gex_only=False
    )

    adata.obs['sample'] = 'B646'
    bdata.obs['sample'] = 'B647'

    fdata = sc.concat([adata, bdata])
    fdata.obs_names_make_unique()

    name_dict = adata.var["feature_types"].to_dict()
    fdata.var["feature_types"] = fdata.var.index.map(name_dict)

    fdata = remove_mito_ribo_hk_lnc_genes(fdata)
    fdata.layers['counts'] = fdata.X.copy()

    sc.pp.filter_cells(fdata, min_genes=200)
    sc.pp.filter_genes(fdata, min_cells=3)

    cleandata = fdata[fdata.obs.log1p_total_counts >= 7, :]
    cleandata = cleandata[cleandata.obs.log1p_n_genes_by_counts >= 6.5, :]

    cleandata = cleandata[cleandata.obs.n_genes_by_counts < 6000, :]
    cleandata = cleandata[cleandata.obs.total_counts < 30000, :]

    sc.pp.normalize_total(cleandata, inplace=True, target_sum=1e4)
    sc.pp.log1p(cleandata)

    def combine_guide_replicates(gdata):
        sgs = gdata.var_names.str.split('_', n=1).str[0]
        sgs_grouped = pd.DataFrame(gdata.X.toarray(), columns=gdata.var_names)
        sgs_grouped = sgs_grouped.groupby(sgs, axis=1).sum()

        cgdata = ad.AnnData(sgs_grouped, obs=gdata.obs, var=pd.DataFrame(index=sgs_grouped.columns))
        return cgdata

    sg_mask = cleandata.var.index.str.startswith('sg')
    sg_counts = pd.DataFrame(cleandata[:, sg_mask].X.toarray())
    singlet_index = (sg_counts.apply(lambda x: x == 0, axis=1).sum(axis=1) == 69)
    doublet_index = (sg_counts.apply(lambda x: x == 0, axis=1).sum(axis=1) == 68)
    doublet_filter_index = np.array(
        ((cleandata[doublet_index][:, -70:].X.sum(axis=1) -
          cleandata[doublet_index][:, -70:].X.max(axis=1)) /
         cleandata[doublet_index][:, -70:].X.sum(axis=1)) < 0.2
    ).flatten()

    map_dict = pd.DataFrame(
        cleandata[:, cleandata.var.feature_types == 'CRISPR Guide Capture'].var_names
    ).to_dict()[0]
    perturb = pd.DataFrame(
        cleandata[:, cleandata.var.feature_types == 'CRISPR Guide Capture'].X.argmax(axis=1)
    )[0].map(map_dict)

    cleandata.obs["perturb"] = "Duplicate"
    cleandata.obs["perturb"][singlet_index.tolist()] = perturb[singlet_index.tolist()].tolist()
    cleandata.obs["perturb"][doublet_index.tolist()][doublet_filter_index] = (
        perturb[doublet_index.tolist()][doublet_filter_index].tolist()
    )

    cdata = combine_guide_replicates(cleandata)
    return cdata

def gex_clustering(cleandata):
    sc.pp.pca(cleandata)
    sc.pp.neighbors(cleandata)
    sc.tl.umap(cleandata)
    sc.tl.leiden(cleandata, resolution=0.2)
    sc.pl.umap(cleandata, color=['n_genes_by_counts', 'leiden'])

    recluster_data = cleandata[(cleandata.obs["leiden"] != '4') & (cleandata.obs["leiden"] != '5')].copy()
    sc.pp.highly_variable_genes(recluster_data, inplace=True, subset=500)
    sc.pp.pca(recluster_data)
    sc.pp.neighbors(recluster_data, n_neighbors=5, n_pcs=0)
    sc.tl.umap(recluster_data)
    return recluster_data

def program_analysis_tf(scdata, tf_linkage, gene_program_names_path="./sc.gene_program.names.csv"):
    """
    Analyze gene programs and assign transcription factor (TF) modules.

    Args:
        scdata: AnnData object containing scRNA-seq data.
        tf_linkage: Hierarchical linkage matrix for TF clustering.
        gene_program_names_path: Path to gene program gene names CSV.

    Returns:
        tf_df (pd.DataFrame): Program scores and module assignments for TF-targeted perturbations.
        chemo_df (pd.DataFrame): Program scores and module assignments for chemokine-targeted perturbations.
        module_means (pd.DataFrame): Program score means for each TF module.
        program_dict (dict): Dictionary of program index to gene lists.
    """

    names_df = pd.read_csv(gene_program_names_path, index_col=0)
    program_dict = {
        '1': names_df['gene'][565:][20:][:-71][:115].tolist(),
        '2': names_df['gene'][565:][20:][:-71][115:176].tolist(),
        '3': names_df['gene'][565:][20:][:-71][176:219].tolist(),
        '4': names_df['gene'][565:][20:][:-71][219:244].tolist(),
        '5': names_df['gene'][565:][20:][:-71][265:].tolist()
    }

    # Build initial DF with all program genes and 'perturb_gene'
    program_genes_all = (
        program_dict['1'] + program_dict['2'] + program_dict['3'] +
        program_dict['4'] + program_dict['5']
    )
    score_mat = pd.DataFrame(
        scdata[:, program_genes_all].X.toarray(),
        index=scdata.obs_names,
        columns=program_genes_all,
    )
    score_mat['perturb_gene'] = scdata.obs['perturb_gene']
    score_df = score_mat[score_mat['perturb_gene'] != 'Duplicate'] \
        .groupby('perturb_gene', as_index=False).mean().dropna().set_index('perturb_gene')

    # Normalize by sgNon-targeting
    score_df = score_df.div(score_df.loc['sgNon-targeting', :], axis=1)
    score_df.dropna(axis=1, inplace=True)

    # Calculate program scores only for genes that exist in the data
    for i in range(5):
        program_genes = [gene for gene in program_dict[str(i+1)] if gene in score_df.columns]
        if program_genes:
            score_df[f'Program {i+1}'] = score_df[program_genes].mean(axis=1)

    # Drop individual gene columns, keep just program and non-NTC perturbs
    score_df = score_df.drop(
        columns=[
            col for col in score_df.columns
            if col in program_genes_all
        ]
    )
    score_df = score_df.loc[score_df.index != 'sgNon-targeting']

    chemokine_list = [
        "sgCcr1", "sgCcr2", "sgCcr4", "sgCcr5", "sgCcr6",
        "sgCcr7", "sgCcr10", "sgCxcr1", "sgCxcr2", "sgCxcr3", "sgCxcr4",
        "sgAckr3", "sgCxcr5", "sgCxcr6", "sgCd74", "sgGpr35"
    ]
    chemo_df = score_df.loc[score_df.index.intersection(chemokine_list)]
    tf_df = score_df.loc[~score_df.index.isin(chemokine_list)]

    # Assign TF modules
    modules = fcluster(tf_linkage, 3, criterion='maxclust')
    tf_df = tf_df.copy()
    tf_df['Module'] = modules

    module_means = tf_df.groupby('Module').mean()

    return tf_df, chemo_df, module_means, program_dict

def module_analysis_chemokine(tf_df, chemo_df, tf_linkage, chem_linkage, scdata, n_modules=3, print_output=True):
    """
    Perform module analysis for transcription factors and chemokines.

    Args:
        tf_df (pd.DataFrame): DataFrame for transcription factors with perturbations as index.
        chemo_df (pd.DataFrame): DataFrame for chemokine perturbations as index.
        tf_linkage: Linkage matrix for TF clustering.
        chem_linkage: Linkage matrix for chemo clustering.
        n_modules (int): Number of clusters/modules to identify.
        print_output (bool): Optionally print the resulting DataFrames.

    Returns:
        module_guides (dict): Maps module number to newline-joined guides in that module (for TF).
        tf_df (pd.DataFrame): tf_df with updated 'Module' column.
        chemo_df (pd.DataFrame): chemo_df with updated 'Module' column.
    """

    # Assign modules to tf_df
    modules_tf = fcluster(tf_linkage, n_modules, criterion='maxclust')
    tf_df = tf_df.copy()
    tf_df['Module'] = modules_tf

    # Generate mapping for modules -> guides (for TF modules)
    module_guides = {}
    for module in tf_df['Module'].unique():
        guides = tf_df[tf_df['Module'] == module].index.tolist()
        module_guides[module] = '\n'.join(guides)

    # Assign modules to chemo_df
    modules_chemo = fcluster(chem_linkage, n_modules, criterion='maxclust')
    chemo_df = chemo_df.copy()
    chemo_df['Module'] = modules_chemo

    if print_output:
        print(tf_df)
        print(chemo_df)

    chemo_gene = chemo_df.index.str.lstrip('sg')
    chemo_gene = chemo_gene[chemo_gene.isin(scdata.var_names)]

    tf_gene = tf_df.index.str.lstrip('sg')
    tf_gene = tf_gene[tf_gene.isin(scdata.var_names)]

    chemo_tf_regdf = pd.concat([pd.DataFrame(scdata[:, chemo_gene].X.toarray(), index=scdata.obs_names, columns=chemo_gene),
                            scdata.obs['perturb_gene']], axis=1)
    chemo_tf_regdf = chemo_tf_regdf[chemo_tf_regdf['perturb_gene'] != 'Duplicate'].groupby('perturb_gene', as_index=False).mean().dropna().set_index('perturb_gene')

    tf_chemo_regdf = pd.concat([pd.DataFrame(scdata[:, tf_gene].X.toarray(), index=scdata.obs_names, columns=tf_gene),
                            scdata.obs['perturb_gene']], axis=1)
    tf_chemo_regdf = tf_chemo_regdf[tf_chemo_regdf['perturb_gene'] != 'Duplicate'].groupby('perturb_gene', as_index=False).mean().dropna().set_index('perturb_gene')
    for i in range(3):
        module_genes = tf_df.groupby('Module').get_group(i+1).index.str.lstrip('sg').tolist()
        module_genes = [gene for gene in module_genes if gene in tf_chemo_regdf.columns]
        tf_chemo_regdf[f'TF Module {i+1}'] = tf_chemo_regdf[module_genes].mean(axis=1)
    tf_chemo_regdf = tf_chemo_regdf.drop(columns=[col for col in tf_chemo_regdf.columns if col in tf_gene])
    for i in range(3):
        module_genes = chemo_df.groupby('Module').get_group(i+1).index.str.lstrip('sg').tolist()
        module_genes = [gene for gene in module_genes if gene in chemo_tf_regdf.columns]
        chemo_tf_regdf[f'Chemokine Module {i+1}'] = chemo_tf_regdf[module_genes].mean(axis=1)
    chemo_tf_regdf = chemo_tf_regdf.drop(columns=[col for col in chemo_tf_regdf.columns if col in chemo_gene])

    chemo_modules = chemo_df['Module'].to_dict()
    tf_modules = tf_df['Module'].to_dict()
    chemo_tf_regdf = pd.DataFrame({
        'Chemokine Module 1': chemo_tf_regdf.loc[[k for k,v in chemo_modules.items() if v==1], ].mean().tolist(),
        'Chemokine Module 2': chemo_tf_regdf.loc[[k for k,v in chemo_modules.items() if v==2], ].mean().tolist(), 
        'Chemokine Module 3': chemo_tf_regdf.loc[[k for k,v in chemo_modules.items() if v==3], ].mean().tolist()
    }, index=['TF Module 1', 'TF Module 2', 'TF Module 3'])

    tf_chemo_regdf = pd.DataFrame({
        'TF Module 1': tf_chemo_regdf.loc[[k for k,v in tf_modules.items() if v==1], ].mean().tolist(),
        'TF Module 2': tf_chemo_regdf.loc[[k for k,v in tf_modules.items() if v==2], ].mean().tolist(),
        'TF Module 3': tf_chemo_regdf.loc[[k for k,v in tf_modules.items() if v==3], ].mean().tolist()
    }, index=['Chemokine Module 1', 'Chemokine Module 2', 'Chemokine Module 3'])

    tf_chemo_regdf = tf_chemo_regdf.melt(ignore_index=False).reset_index()
    chemo_tf_regdf = chemo_tf_regdf.melt(ignore_index=False).reset_index()

    tf_chemo_regdf.columns = ['Chemokine Module', 'TF Module', 'Score']
    chemo_tf_regdf.columns = ['TF Module', 'Chemokine Module', 'Score']

    return tf_chemo_regdf, chemo_tf_regdf