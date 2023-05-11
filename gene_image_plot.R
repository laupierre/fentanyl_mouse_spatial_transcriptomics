library (Seurat)
library (ggpubr)


### Slide 1
brain <- readRDS ("brain_slide1_G1G3_groups.rds")

Images (brain)
# G1-1A, G3-1B, G1-1C, G3-1D
# "slice1"   "slice1.2" "slice1.3" "slice1.4"

p1c <- SpatialDimPlot(brain, images=c("slice1", "slice1.3"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")
p2c <- SpatialDimPlot(brain, images=c("slice1.2", "slice1.4"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")

p <- ggarrange(p1c, p2c, nrow=2, labels=c("slices G1", "slices G3"))
p
ggsave ("slide1.pdf")



### Slide 2
brain2 <- readRDS ("brain_slide2_G2G4_groups.rds")

Images (brain2)
# G2-2A, G4-2B, G2-2C, G4-2D
# "slice1"   "slice1.2" "slice1.3" "slice1.4"

p3c <- SpatialDimPlot(brain2, images=c("slice1", "slice1.3"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")
p4c <- SpatialDimPlot(brain2, images=c("slice1.2", "slice1.4"), label=FALSE, label.size= 5, repel = TRUE) + theme(legend.position = "none")

p <- ggarrange(p3c, p4c, nrow=2, labels=c("slices G2", "slices G4"))
p
ggsave ("slide2.pdf")


## background
pa1 <- ggarrange(p1c, p2c, nrow=2, labels=c("G1", "G3"))
pa2 <- ggarrange(p3c, p4c, nrow=2, labels=c("G2", "G4"))
pa3 <- pa1 | pa2
pa3
ggsave ("background plot.pdf", pa3)



## gene
gene <- "Ttr"

# G1 group
p1 <- SpatialFeaturePlot(brain, images=c("slice1", "slice1.3"), features = gene)
# G3 group
p2 <- SpatialFeaturePlot(brain, images=c("slice1.2", "slice1.4"), features = gene)


# alpha = c(0.3, 3)
# G2 group
p3 <- SpatialFeaturePlot(brain2, images=c("slice1", "slice1.3"), features = gene)
# G4 group
p4 <- SpatialFeaturePlot(brain2, images=c("slice1.2", "slice1.4"), features = gene)


pa1 <- ggarrange(p1, p2, nrow=2, labels=c("G1", "G3"))
pa2 <- ggarrange(p3, p4, nrow=2, labels=c("G2", "G4"))
pa3 <- pa1 | pa2
pa3
ggsave (paste (gene, " plot.pdf", sep=""), pa3)





