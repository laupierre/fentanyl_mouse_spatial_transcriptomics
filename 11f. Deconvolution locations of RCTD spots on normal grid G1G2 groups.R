library (Seurat)
library (ggpubr)
library (sctransform)
library (RColorBrewer)

### This is for 4 slices together (G1 and G2 groups)

### Pseudoslide
brain <- readRDS ("brain_G2G1_groups.rds")

table (gsub ("_.*", "", row.names (brain@meta.data)))
#  1A   1C   2A   2C 
#3276 3552 3336 3361 



meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G2-2C/cell_proportions_allen.txt")
colnames (meta)[1] <- "Barcode"
meta$Barcode <- paste ("2C_", meta$Barcode, sep="")
colnames (meta)[2] <- "Hip"
meta <- meta[meta$Hip != "", ]

meta <- meta[meta$Barcode %in%  row.names (brain@meta.data), ]
idx <- match (meta$Barcode, row.names (brain@meta.data))
brain@meta.data$location[idx] <- meta$Hip


rm (meta)
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G2-2A/cell_proportions_allen.txt")
colnames (meta)[1] <- "Barcode"
meta$Barcode <- paste ("2A_", meta$Barcode, sep="")
colnames (meta)[2] <- "Hip"
meta <- meta[meta$Hip != "", ]

meta <- meta[meta$Barcode %in%  row.names (brain@meta.data), ]
idx <- match (meta$Barcode, row.names (brain@meta.data))
brain@meta.data$location[idx] <- meta$Hip


rm (meta)
meta <-  read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1C/cell_proportions_allen.txt")
colnames (meta)[1] <- "Barcode"
meta$Barcode <- paste ("1C_", meta$Barcode, sep="")
colnames (meta)[2] <- "Hip"
meta <- meta[meta$Hip != "", ]

meta <- meta[meta$Barcode %in%  row.names (brain@meta.data), ]
idx <- match (meta$Barcode, row.names (brain@meta.data))
brain@meta.data$location[idx] <- meta$Hip


rm (meta)
meta <-  read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A/cell_proportions_allen.txt")
colnames (meta)[1] <- "Barcode"
meta$Barcode <- paste ("1A_", meta$Barcode, sep="")
colnames (meta)[2] <- "Hip"
meta <- meta[meta$Hip != "", ]

meta <- meta[meta$Barcode %in%  row.names (brain@meta.data), ]
idx <- match (meta$Barcode, row.names (brain@meta.data))
brain@meta.data$location[idx] <- meta$Hip



brain@meta.data$location <- gsub ("=.*", "", brain@meta.data$location)
brain@meta.data$allen <- brain@meta.data$location
table (brain@meta.data$allen)


keep <- c("CA1-ProS", "CA2-IG-FC", "CA3", "DG")
brain@meta.data$allen [which (!brain@meta.data$allen %in% keep)] <- "unknown"


colors= c("CA1-ProS" = "red", "CA2-IG-FC" = "pink", "CA3" = "green", "DG" = "blue", "unknown" = "grey")

p1 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1") + ggplot2::theme(legend.position = "none")
p2 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1.2") + ggplot2::theme(legend.position = "none")
p3 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1.3") + ggplot2::theme(legend.position = "none")
p4 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1.4") + ggplot2::theme(legend.position = "none")

pa1 <- ggarrange (p2 | p1, nrow=1, labels="G2")
pa2 <- ggarrange (p4 | p3, nrow=1, labels="G1")
pall1 <- ggarrange (pa1, pa2, nrow=2)
pall1 

ggsave ("Deconvolution locations of RCTD spots substructures plot.jpeg", pall1, units="in", width=8, height=8, dpi=300)




