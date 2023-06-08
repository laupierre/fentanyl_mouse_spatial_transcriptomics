library (Seurat)
library (ggpubr)
library (sctransform)
library (RColorBrewer)


### Pseudoslide
brain <- readRDS ("brain_G2G1_groups.rds")

table (gsub ("_.*", "", row.names (brain@meta.data)))
#  1A   1C   2A   2C 
#3276 3552 3336 3361 


meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G2_2C_sniv02.csv", sep=",")
meta$Barcode <- paste ("2C_", meta$Barcode, sep="")
colnames (meta)[1] <- "Barcode"
colnames (meta)[3] <- "Hip"
meta <- meta[meta$Hip != "", ]

idx <- match (meta$Barcode, row.names (brain@meta.data))
stopifnot (length (idx) != 0)
brain@meta.data$location[idx] <- meta$Hip


rm (meta)
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G2_2A_sniv01.csv", sep=",")
meta$Barcode <- paste ("2A_", meta$Barcode, sep="")
colnames (meta)[1] <- "Barcode"
colnames (meta)[3] <- "Hip"
meta <- meta[meta$Hip != "", ]

idx <- match (meta$Barcode, row.names (brain@meta.data))
stopifnot (length (idx) != 0)
brain@meta.data$location[idx] <- meta$Hip


rm (meta)
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1C_shmv02.csv", sep=",")
meta$Barcode <- paste ("1C_", meta$Barcode, sep="")
colnames (meta)[1] <- "Barcode"
colnames (meta)[3] <- "Hip"
meta <- meta[meta$Hip != "", ]

idx <- match (meta$Barcode, row.names (brain@meta.data))
stopifnot (length (idx) != 0)
brain@meta.data$location[idx] <- meta$Hip


rm (meta)
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01.csv", sep=",")
meta$Barcode <- paste ("1A_", meta$Barcode, sep="")
colnames (meta)[1] <- "Barcode"
colnames (meta)[3] <- "Hip"
meta <- meta[meta$Hip != "", ]

idx <- match (meta$Barcode, row.names (brain@meta.data))
stopifnot (length (idx) != 0)
brain@meta.data$location[idx] <- meta$Hip


Images (brain)

colors= c(CA1 = "red", CA2 = "pink", "CA3" = "green", DG= "blue", unknown = "grey")

p1 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1") + ggplot2::theme(legend.position = "none")
p2 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1.2") + ggplot2::theme(legend.position = "none")
p3 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1.3") + ggplot2::theme(legend.position = "none")
p4 <- SpatialDimPlot(brain, cols=colors, group.by="location", images = "slice1.4") + ggplot2::theme(legend.position = "none")

pa1 <- ggarrange (p2 | p1, nrow=1, labels="G2")
pa2 <- ggarrange (p4 | p3, nrow=1, labels="G1")
pall1 <- ggarrange (pa1, pa2, nrow=2)

ggsave ("selected cells highlight substructures plot.pdf", pall1, width=8, height=8)
# ggsave ("selected cells highlight substructures plot.jpeg", p, units="in", width=8, height=8, dpi=300)





