library (Seurat)
library (ggpubr)


### Slide 1
brain <- readRDS ("brain_slide1_G1G3_groups.rds")

Images (brain)
# G1-1A, G3-1B, G1-1C, G3-1D
# "slice1"   "slice1.2" "slice1.3" "slice1.4"

p1c <- SpatialDimPlot(brain, images=c("slice1", "slice1.3"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")
p2c <- SpatialDimPlot(brain, images=c("slice1.2", "slice1.4"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")

ggarrange(p1c, p2c, nrow=2, labels=c("slices G1", "slices G3"))





### Slide 2
brain2 <- readRDS ("brain_slide2_G2G4_groups.rds")

Images (brain2)
# G2-2A, G4-2B, G2-2C, G4-2D
# "slice1"   "slice1.2" "slice1.3" "slice1.4"

p1c <- SpatialDimPlot(brain2, images=c("slice1", "slice1.3"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")
p2c <- SpatialDimPlot(brain2, images=c("slice1.2", "slice1.4"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")

ggarrange(p1c, p2c, nrow=2, labels=c("slices G2", "slices G4"))



gene <- "Arc"

# alpha = c(0.3, 3)
# G2 group
p3 <- SpatialFeaturePlot(brain2, images=c("slice1", "slice1.3"), features = gene)
# G4 group
p3 <- SpatialFeaturePlot(brain2, images=c("slice1.2", "slice1.4"), features = gene)


ggarrange(p1c, p2c, p3c, p4c, nrow=2, ncol=2, labels=c("G2", "G4"))

# p1 | p2 | p3 | p4
ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, labels=c("G2", "G4"))






