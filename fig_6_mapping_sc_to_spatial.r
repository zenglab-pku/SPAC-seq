library(tidyverse)

options(stringsAsFactors = F)
library(tidyverse)
library(CellTrek)
library(viridis)
library(ConsensusClusterPlus)

# Read in the sc data
sc_meta <- read.table("sc.reclustered.obs.csv", sep=",", header=T, row.names=1)
cdata <- CreateSeuratObject(counts=t(sc_matrix), meta.data=sc_meta)

cdata <- FindVariableFeatures(cdata)
cdata$group <- cdata@meta.data$perturb_gene
Idents(cdata) <- 'group'

# mt_genes is a list of t signature genes
cv_genes <- mt_genes[mt_genes %in% rownames(cdata)]
t_high_data <- subset(cdata, features = cv_genes)

t_high_data <- NormalizeData(t_high_data)
t_high_data <- ScaleData(t_high_data)
t_high_data <- RunPCA(t_high_data)
t_high_data <- RunUMAP(t_high_data, reduction='pca', dims=1:20)

DimPlot(t_high_data, label=T, group.by="perturb_gene", label.size = 4.5)

# Read in the spatial data
spatial_meta <- read.table("spatial.reclustered.obs.csv", sep=",", header=T)
rownames(spatial_matrix) <- rownames(spatial_meta)
spatial_data <- CreateSeuratObject(counts=t(spatial_matrix), meta.data=spatial_meta)

spatial_coords <- read.table("spatial.reclustered.spatial.txt", sep="\t")
spatial_data@meta.data$coord_x <- spatial_coords[, 1]
spatial_data@meta.data$coord_y <- spatial_coords[, 2]

spatial_data@images$Spatial <- new(
    Class="VisiumV1",
    assay="RNA",
    key="image_",
    coordinates=data.frame(
        imagerow = spatial_data@meta.data$coord_x,
        imagecol = spatial_data@meta.data$coord_y
    )%>%magrittr::set_rownames(rownames(spatial_data@meta.data)),
    scale.factors=scalefactors(spot = 138.656, fiducial = 223.9828, hires = 0.1139861, lowres = 0.03419583)
)

spatial_data <- RenameCells(spatial_data, new.names = make.names(Cells(spatial_data)))
cdata <- RenameCells(cdata, new.names=make.names(Cells(cdata)))

# sample 1
tdata <- subset(spatial_data, subset = marker == 'A')
cdata@meta.data$type <- 'T'
tdata <- NormalizeData(tdata)
traint <- CellTrek::traint(st_data = tdata, sc_data = cdata, sc_assay = "RNA", st_assay = "RNA", cell_names = 'type', nfeatures = 2000)
DimPlot(traint, group.by='type', raster=F)

celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=cdata, sc_assay = 'RNA',
                                intp_pnt=500, intp_lin=F, nPCs=10, ntree=200, 
                                dist_thresh=0.55, top_spot=4, spot_n=10,
                                reduction='pca', intp=T, repel_r=20, repel_iter=20, keep_model=T)$celltrek

# custom_colors is a list of colors for each guide
generate_colors <- colorRampPalette(custom_colors)(35)
  
options(repr.plot.width = 10, repr.plot.height = 5)
ggplot(celltrek@meta.data %>% filter(perturb_gene != 'Duplicate'), aes(x = coord_x, y = coord_y, color = .data[['perturb_gene']])) +
geom_point(size=0.5) +
labs(title = paste("Spatial Distribution of perturb"),
      x = "X Coordinate",
      y = "Y Coordinate",
      color = 'perturb_gene') +
theme_minimal()+
scale_color_manual(values=generate_colors)

c_tab <- celltrek@meta.data %>% dplyr::select(c('coord_x', 'coord_y', 'perturb', 'perturb_gene', 'phenotype'))
write.table(c_tab, 'sc.double.t_map.csv', sep='\t', quote=FALSE)