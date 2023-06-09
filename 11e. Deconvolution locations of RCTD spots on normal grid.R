library (Seurat)
library (ggplot2)

options(Seurat.object.assay.version = "v5")

## This is only for 1 slice

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A/cell_proportions_allen.txt")


sample.name <- gsub (".*/", "", data.dir)
brain <- Load10X_Spatial (data.dir, filename= "filtered_feature_bc_matrix.h5")
brain[["Spatial"]] <- as(brain[["Spatial"]], Class = "Assay5")

# Remove cells with mitochondrial contamination
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "percent_mito")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")
brain <- brain[ ,brain$nFeature_Spatial > 500 & brain$percent_mito < 20 & brain$percent_hb < 20]

# remove mitochondrial genes
brain <- brain[!grepl("^mt-", rownames(brain)), ]
# remove Hemoglobin genes (if that is a problem on your data)
brain <- brain[!grepl("^Hb.*-", rownames(brain)), ]

colnames (meta)[1] <- "Barcode"
meta1 <- meta

tokeep <- intersect (meta1$Barcode, row.names (brain@meta.data))
brain <- subset (brain, cells = tokeep)
dim (brain[["Spatial"]]$counts)

meta1 <- meta1[meta1$Barcode %in% tokeep, ]
idx <- match (tokeep, row.names (brain@meta.data))
#brain@meta.data$allen <- meta1$Barcode[idx]
#table (brain@meta.data$allen == row.names (brain@meta.data))
brain@meta.data$allen <- meta1$cell_type1[idx]
brain@meta.data$allen <- gsub ("=.*","", brain@meta.data$allen)
brain@meta.data$group <- sample.name
brain@meta.data$cell <- paste (row.names(brain@meta.data), sample.name, sep="-")

brain <- UpdateSeuratObject(brain)

Idents (brain) <- "allen"
p1 <- SpatialDimPlot(brain, label = TRUE, label.size = 3, stroke=NA) + theme(legend.position='none')
p1
#SpatialDimPlot(brain, group.by = c("allen"))


## subset to hippocampus (this mess up with the picture resolution)
#hip <- subset(brain, idents = c("CA1-ProS", "CA2-IG-FC", "CA3", "DG"))
#p2 <- SpatialDimPlot(hip, crop=FALSE, label = FALSE, pt.size.factor = 0.8, alpha = 0.9) + ggplot2::theme(legend.position = "none")
##ggsave("preselected_cells.jpg", height=8, width=8, units='in', dpi=300)
#ggplot2::ggsave("preselected_cells_normal_slice.pdf", height=8, width=8)

keep <- c("CA1-ProS", "CA2-IG-FC", "CA3", "DG")
brain@meta.data$allen [which (!brain@meta.data$allen %in% keep)] <- "unknown"
colors= c("CA1-ProS" = "blue", "CA2-IG-FC" = "red", "CA3" = "pink", "DG"= "darkgreen", "unknown" = "white")

p2 <- SpatialDimPlot(brain, cols=colors, group.by="allen", images = "slice1") + ggplot2::theme(legend.position = "none")
p2
ggsave ("version 1 plot.pdf", p2, width=8, height=8)


## Add original annotation
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01.csv", sep=",")
colnames (meta)[2] <- "Barcode"
colnames (meta)[3] <- "Hip"
meta <- meta[meta$Hip != "", ]

idx <- match (meta$Barcode, row.names (brain@meta.data))
stopifnot (length (idx) != 0)
brain@meta.data$preselection <- "Unknown"
brain@meta.data$preselection[idx] <- meta$Hip

colors= c(CA1 = "red", CA2 = "pink", "CA3" = "darkgreen", DG= "blue", unknown = "grey")

p3 <- SpatialDimPlot(brain, cols=colors, group.by="preselection", images = "slice1") + ggplot2::theme(legend.position = "none")
p3
ggsave ("version 2 plot.pdf", p3, width=8, height=8)






######################
## get the palo colors
## https://winnie09.github.io/Wenpin_Hou/pages/Palo.html

#saveRDS (brain, "brain_palo.rds")
#brain <- readRDS ("brain_palo.rds")

library(Palo)

# Get the spatial coordinates (lowres is by default) (or get the UMAP coordinates when available, ie: d[[“umap”]]@cell.embeddings)
coords <- brain@images$slice1@coordinates
coords <- coords[ ,c("row", "col")]

# Get the spot clusters
Idents (brain) <- "allen"
#cl <- Idents(brain)
cl <- Idents(brain)
names (cl) <- row.names (coords)
cl <- factor (cl)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

pal <- gg_color_hue(length(unique(cl)))

## FIXME
# palopal <- Palo (coords, cl, pal)
palopal <- pal

SpatialDimPlot(brain, label = TRUE, label.size = 3, stroke=NA) + theme(legend.position='none') +
               scale_fill_manual(values=palopal)






