library(fmsb)
library(RColorBrewer)
library(scales)

df <- read.csv("./cibersort/results.csv", row.names=1)
df <- na.omit(df)
zscore_df <- as.data.frame(scale(df))
colnames(zscore_df)
cell_types <- zscore_df[,c('Plasma.cells',"T.cells.CD8",'T.cells.follicular.helper','T.cells.regulatory..Tregs.','NK.cells.activated',
        "Monocytes",'Macrophages.M1','Macrophages.M2','Dendritic.cells.activated')]
ROI_df <- cell_types[c("sgnon-targeting","sgCd276","sgIcam1"),]
ROI_df
ROI_df <- rbind(
  max = rep(2, ncol(ROI_df)),
  min = rep(-2, ncol(ROI_df)),
  ROI_df
)
options(repr.plot.width = 10, repr.plot.height = 9)

radarchart(
  ROI_df,
  axistype = 1,
  pcol = c("#00bfff", "#e68ab8", "#f0e68c"),
  caxislabels = c("-2", "-1", "0", "1", "2"),
  plwd = 3,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 1.5,
  axislabcol = "black",
  vlabels = colnames(data),
  vlcex = 1.2,
  calcex = 1.2
)

legend(x = "bottom", legend = c("sgnon-targeting","sgCd276","sgIcam1"), horiz = T,
       bty = "n", pch = 16, col = c("#00bfff", "#e68ab8", "#f0e68c"), text.col = "black",
       cex = 1, pt.cex = 2)

pdf("ROIs_cell_type_proportions.pdf", width = 10, height = 9)

radarchart(
  ROI_df,
  axistype = 1,
  pcol = c("#00bfff", "#e68ab8", "#f0e68c"),
  caxislabels = c("-2", "-1", "0", "1", "2"),
  plwd = 3,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 1.5,
  axislabcol = "black",
  vlabels = colnames(data),
  vlcex = 1.2,
  calcex = 1.2
)

legend(x = "bottom", legend = c("sgnon-targeting","sgCd276","sgIcam1"), horiz = T,
       bty = "n", pch = 16, col = c("#00bfff", "#e68ab8", "#f0e68c"), text.col = "black",
       cex = 1, pt.cex = 2)
       
dev.off()
ROI_df <- cell_types[c("sgCcn1","sgCd99"),]
ROI_df <- rbind(
  max = rep(2, ncol(ROI_df)),
  min = rep(-2, ncol(ROI_df)),
  ROI_df
)
options(repr.plot.width = 10, repr.plot.height = 9)

radarchart(
  ROI_df,
  axistype = 1,
  pcol = c("#ffb366","#e32636"),
  pfcol = alpha(c("#ffb366","#e32636"), 0.6),
  caxislabels = c("-2", "-1", "0", "1", "2"),
  plwd = 3,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 1.5,
  axislabcol = "black",
  vlabels = colnames(data),
  vlcex = 1.2,
  calcex = 1.2
)

legend(x = "bottom", legend = c("sgCcn1","sgCd99"), horiz = T,
       bty = "n", pch = 16, col = c("#ffb366","#e32636"), text.col = "black",
       cex = 1, pt.cex = 2)

pdf("sgCcn1_cell_type_proportions.pdf", width = 10, height = 9)

radarchart(ROI_df, axistype = 1, pcol = c("#ffb366"), pfcol = alpha(c("#ffb366"), 0.6),
      caxislabels = c("-2", "-1", "0", "1", "2"), plwd = 3, plty = 1, cglcol = "grey",
      cglty = 1, cglwd = 1.5, axislabcol = "black", vlabels = colnames(data),
      vlcex = 1.2, calcex = 1.2
)
       
dev.off()
ROI_df <- cell_types[c("sgCd99"),]
ROI_df <- rbind(
  max = rep(2, ncol(ROI_df)),
  min = rep(-2, ncol(ROI_df)),
  ROI_df
)

options(repr.plot.width = 10, repr.plot.height = 9)

radarchart(
  ROI_df,
  axistype = 1,
  pcol = c("#e32636"),
  pfcol = alpha(c("#e32636"), 0.6),
  caxislabels = c("-2", "-1", "0", "1", "2"),
  plwd = 3,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 1.5,
  axislabcol = "black",
  vlabels = colnames(data),
  vlcex = 1.2,
  calcex = 1.2
)

pdf("sgCd99_cell_type_proportions.pdf", width = 10, height = 9)

radarchart(
  ROI_df,
  axistype = 1,
  pcol = c("#e32636"),
  pfcol = alpha(c("#e32636"), 0.6),
  caxislabels = c("-2", "-1", "0", "1", "2"),
  plwd = 3,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 1.5,
  axislabcol = "black",
  vlabels = colnames(data),
  vlcex = 1.2,
  calcex = 1.2
)
       
dev.off()