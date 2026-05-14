library(tidyverse)
library(ggrepel)

deseq2_out <- read.csv('./bulkdata/Group_CD44PBS_vs_interPBS_DESeq2.txt', sep = '\t')
pbs_results <- deseq2_out

deseq2_out_sorted <- deseq2_out %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

tab <- deseq2_out_sorted

tab$padj <- ifelse(tab$padj < 1e-50, 1e-50, tab$padj)
tab$log2FoldChange <- ifelse(tab$log2FoldChange > 3.5, 3.5, tab$log2FoldChange)
tab$log2FoldChange <- ifelse(tab$log2FoldChange < -3.5, -3.5, tab$log2FoldChange)

# Assign colors based on fold change and significance
tab$color <- 'Not Significant'
tab$color[tab$log2FoldChange > 0.5 & tab$padj < 0.001] <- 'Up'
tab$color[tab$log2FoldChange < -0.5 & tab$padj < 0.001] <- 'Down'
tab$color <- factor(tab$color, levels = c('Down', 'Not Significant', 'Up'))

# Create the plot with soft colors and labeled genes
p <- ggplot(tab, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color), size = 1.5) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  theme_classic() +
  labs(
    title = "sgCd44 ~ sgIntergenic PBS",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = 'Significance'
  ) +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(1, 50)) +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray60") +
  geom_label_repel(
    data = subset(tab, (padj < 0.01) & (abs(log2FoldChange) > 0.5)),
    aes(label = Gene_name),
    size = 3, fill = "white",
    box.padding = 0.4,
    label.size = 0,
    label.r = unit(0.2, "lines"), # Rounded corners
    color = "black"
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.3),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.line = element_line(size = 0.5, color = "gray40"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

print(p)

gene_list <- tab %>% filter(padj < 0.05) %>% pull(log2FoldChange)
names(gene_list) <- tab %>% filter(padj < 0.05) %>% pull(Gene_name)
head(gene_list)

library(msigdbr)

# Get MSigDB immune-related gene sets for mouse (C7 category)
immune_sets <- msigdbr(species = "Mus musculus", category = "C7") %>%
  dplyr::select(gs_name, gene_symbol)  # Ensure you're using gene symbols

# Prepare geneList and ensure it is sorted
geneList <- sort(gene_list, decreasing = TRUE)

# Run GSEA
gsea_results <- GSEA(geneList, TERM2GENE = immune_sets, pvalueCutoff = 2, minGSSize = 1, maxGSSize = 1000)

id = gsea_results@result$ID == 'GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_UP'
# gseaplot2(gsea_results, geneSetID = gsea_results@result$ID[id], title = gsea_results@result$Description[id])

p1 <- gseaplot2(
  gsea_results,
  geneSetID = gsea_results@result$ID[id],
  title = paste('GOLDRATH_STEM-LIKE_VS_MEMORY_CD8_TCELL_UP', "\nNES:", round(gsea_results@result$NES[id], 2),
    "\nFDR:", gsea_results@result$qvalue[id]),
  color = "green"
)

print(p1)