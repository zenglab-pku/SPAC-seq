library(nichenetr)
library(tidyverse)

lr_network <- readRDS('./RNA/nichenet/lr_network_mouse_21122021.rds')
ligand_target_matrix <- readRDS('./RNA/nichenet/ligand_target_matrix_nsga2r_final_mouse.rds')
weight_network <- readRDS('./RNA/nichenet/weighted_networks_nsga2r_final_mouse.rds')

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
df <- read_csv('./cluster_specfic_cell_type_expression.csv', col_names = TRUE)
head(df)
cdf <- df
cdf$cluster <- rep(0:11, 14)

cdf$sample <- paste(cdf$annotations, cdf$cluster, sep = "_")
row.names(cdf) <- cdf$sample

cdf <- t(cdf)
row.names(cdf) <- convert_alias_to_symbols(row.names(cdf), 'mouse', verbose=TRUE)

cdf <- cdf[1:(nrow(cdf) - 3), ]

sinfo <- data.frame(sample = colnames(cdf))
sinfo$cluster <- sub(".*_", "", sinfo$sample)
receptors <- lr_network %>% pull(to) %>% unique()
expressed_genes_receiver <- t(cdf[, grepl("T_[0-9]*$", colnames(cdf))]) %>%
    apply(2, function(x){10 * (2^as.numeric(x) - 1)}) %>%
    apply(2, function(x){log2(mean(x) + 1)}) %>%
    .[. >= 500] %>% names()
expressed_receptors <- intersect(receptors, expressed_genes_receiver)
print(c(length(expressed_receptors), 'Cd44' %in% expressed_receptors))
ligands <- lr_network %>% pull(from) %>% unique()
expressed_genes_sender <- t(cdf[, grepl("Macrophages_[0-9]*$", colnames(cdf))]) %>%
    apply(2, function(x){10 * (2^as.numeric(x) - 1)}) %>%
    apply(2, function(x){log2(mean(x) + 1)}) %>%
    .[. >= 500] %>% names()
expressed_ligands <- intersect(ligands, expressed_genes_sender)
print(c(length(expressed_ligands), 'Spp1' %in% expressed_ligands))
# Only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
T_marker <- c('Cd8a', 'Cd8b', 'Il7r', 'Cd27', 'Ccr4', 'Eomes', 'Gpr183', 'Cd69', 'Ifng', 'Fgfbp2', 'Gzmh', 'Gzmk', 'Gzma', 'Gzmb', 'Prf1', 'Gnly', 'Nkg7', 'Lag3', 'Pdcd1', 'Tigit', 'Ctla4', 'Tox', 'Sirpg', 'Tnfrsf9', 'Tnfrsf18', 'Cxcr6', 'Cxcl13', 'Tcf7', 'Foxo1')
geneset_oi <- T_marker %>% .[. %in% rownames(ligand_target_matrix)]
length(geneset_oi)
potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
length(potential_ligands)
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
head(ligand_activities)
(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected))))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

best_upstream_ligands

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>% bind_rows()

# Check for missing values and handle them
if (any(is.na(active_ligand_target_links_df$weight))) {
  active_ligand_target_links_df <- active_ligand_target_links_df %>%
    filter(!is.na(weight))  # Remove missing values
}

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Prioritized CAF-ligands", "p-EMT genes in malignant cells",
                    color = "orange", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0) +
  theme(plot.width = unit(20, "in"), plot.height = unit(6, "in"))

p_ligand_target_network
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weight_network$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both")

# Convert vis_ligand_receptor_network to a data frame for select and filter operations
vis_ligand_receptor_network_df <- as.data.frame(vis_ligand_receptor_network)

vis_ligand_receptor_network_df <- vis_ligand_receptor_network_df %>%
  select(-contains("Wnt")) %>%
  filter(!rownames(vis_ligand_receptor_network_df) %in% grep("F", rownames(vis_ligand_receptor_network_df), value = TRUE)) 

ggplot(data = as.data.frame(as.table(as.matrix(vis_ligand_receptor_network_df))), aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "mediumvioletred") +
  labs(y = "Prioritized Ligands", x = "Receptors expressed", fill = "Prior interaction potential") +
  theme(plot.width = unit(9, "in"), plot.height = unit(9, "in"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size=8))  # Adjust font size and rotate xticks

ggsave('./plots/cluster_lr_ligand_receptor.pdf', width = 12, height = 9, device = 'pdf', dpi = 300)

library(reshape2)
melted_vis_ligand_receptor_network_df <- melt(vis_ligand_receptor_network_df)
melted_vis_ligand_receptor_network_df$source <- rownames(vis_ligand_receptor_network_df)
colnames(melted_vis_ligand_receptor_network_df) <- c("target", "weight", "source")

melted_vis_ligand_receptor_network_df <- melted_vis_ligand_receptor_network_df[c(3, 1, 2)]
print(melted_vis_ligand_receptor_network_df)
write.csv(melted_vis_ligand_receptor_network_df, './cluster_lr_ligand_receptor_network.csv', row.names = FALSE)

# Use flows data to draw a chord diagram
# Adjust the proportion of source in flows so that the proportion of source and target nodes is half and half
source_proportion <- table(melted_vis_ligand_receptor_network_df$source) / sum(table(melted_vis_ligand_receptor_network_df$source))
target_proportion <- table(melted_vis_ligand_receptor_network_df$target) / sum(table(melted_vis_ligand_receptor_network_df$target))

# Make sure the names of source_proportion and target_proportion are consistent
all_levels <- union(names(source_proportion), names(target_proportion))
source_proportion <- source_proportion[all_levels]
target_proportion <- target_proportion[all_levels]
source_proportion[is.na(source_proportion)] <- 0
target_proportion[is.na(target_proportion)] <- 0

combined_proportion <- (source_proportion + target_proportion) / 2

# Make sure the length of grid_colors matches the number of unique values in flows$source
grid_colors <- RColorBrewer::brewer.pal(n = min(length(unique(melted_vis_ligand_receptor_network_df$source)), 11), name = "RdBu")
# Use colorRamp function to define edge colors
# Make sure the length of breaks matches that of colors
breaks <- seq(min(melted_vis_ligand_receptor_network_df$weight), max(melted_vis_ligand_receptor_network_df$weight), length.out = 11)
edge_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(length(breaks) - 1)

library(circlize)

# Draw a chord diagram with arrows for edges and increased spacing between source and target
setwd('~/stereoseq/20240502-SPACseq/')
pdf(file = "./chord_diagram.pdf", width = 8, height = 8)
circlize::chordDiagram(melted_vis_ligand_receptor_network_df, transparency = 0.7, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1),
             col = edge_colors[findInterval(melted_vis_ligand_receptor_network_df$weight, breaks)], directional = 1, direction.type = "arrows", link.arr.type = "big.arrow")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(5), sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
dev.off()