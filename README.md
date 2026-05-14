# SPAC-seq analysis scripts

Analysis code for SPAC-seq paper.

* **fig_2_cibersort.r**: Plots z-scored CIBERSORT immune cell fractions (radar charts and PDF export) for selected sgRNAs versus non-targeting controls.

* **fig_3_cluster_guide_analysis.py**: Builds guide-by-cluster matrices, runs chi-square / enrichment-style statistics, and summarizes guide usage per spatial cluster from AnnData objects.

* **fig_3_load_stereoseq.py**: Loads BGI Stereo-seq GEM/bin tables and optional tissue images into `AnnData` with spatial coordinates.

* **fig_3_module_program_analysis.py**: Fits NMF on highly variable genes, runs consensus clustering on NMF components, and derives co-expression / program-level summaries.

* **fig_3_perturbation_preprocessing.py**: Concatenates spatial samples, separates RNA from guide features, merges non-targeting guide channels, and prepares normalized perturbation AnnData for downstream analysis.

* **fig_3_transcriptomic_preprocessing.py**: End-to-end spatial transcriptome workflow (neighbors, clustering, optional scVI/cellcharter-related steps, DE helpers) on binned Stereo-seq objects.

* **fig_3_utils.py**: Shared utilities—gene symbol checks, filtering of mitochondrial/ribosomal/housekeeping genes, batch loading of `.h5ad` files, and common plotting helpers.

* **fig_5_bulk_rna.r**: Reads DESeq2 bulk RNA-seq results and draws a labeled volcano plot for group comparisons.

* **fig_5_delauney_interaction.py**: Uses Squidpy to build spatial neighbor graphs and a binary neighborhood-enrichment matrix from interaction z-scores for categorical annotations.

* **fig_5_interaction_beta_regression.py**: Scores cell-type marker panels on spatial bins and fits regularized (elastic net) models relating RNA features to matched guide capture in a chosen cluster.

* **fig_5_ligand_receptor.r**: Runs NicheNet-style ligand–receptor inference on cluster-resolved cell-type expression (sender vs receiver definitions and activity scores).

* **fig_6_interaction_kde.py**: KDE-based spatial analysis comparing two guides along a gene-expression gradient, with smoothing and a simple group comparison test.

* **fig_6_mapping_sc_to_spatial.r**: Seurat preprocessing of perturb-seq and spatial objects, then CellTrek-style mapping of single cells onto tissue coordinates.

* **fig_6_perturb-seq_analysis.py**: Loads 10x multiome-style matrices with CRISPR features, QC filtering, normalization, clustering, and perturbation-aware single-cell summaries.
