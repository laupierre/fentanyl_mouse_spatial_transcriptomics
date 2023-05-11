library (Seurat)
library (ggpubr)


# Slide 1
brain <- readRDS ("brain_slide1_G1G3_groups.rds")

Images (brain)
# G1-1A, G3-1B, G1-1C, G3-1D
# "slice1"   "slice1.2" "slice1.3" "slice1.4"

p1c <- SpatialDimPlot(brain, images=c("slice1", "slice1.3"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")
p2c <- SpatialDimPlot(brain, images=c("slice1.2", "slice1.4"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")

ggarrange(p1c, p2c, nrow=2, labels=c("slices G1", "slices G3"))





# Slide 2
brain <- readRDS ("brain_slide2_G2G4_groups.rds")

Images (brain)
# G2-2A, G4-2B, G2-2C, G4-2D
# "slice1"   "slice1.2" "slice1.3" "slice1.4"

p1c <- SpatialDimPlot(brain, images=c("slice1", "slice1.3"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")
p2c <- SpatialDimPlot(brain, images=c("slice1.2", "slice1.4"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")

ggarrange(p1c, p2c, nrow=2, labels=c("slices G2", "slices G4"))


gene <- "Arc"

# G2 group
p1 <- SpatialFeaturePlot(brain, images=c("slice1", "slice1.3"), features = gene)
# G4 group
p2 <- SpatialFeaturePlot(brain, images=c("slice1.2", "slice1.4"), features = gene)

# p1 | p2
ggarrange(p1, p2, nrow=2, labels=c("G2", "G4"))






