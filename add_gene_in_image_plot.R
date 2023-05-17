library (Seurat)
library (ggpubr)
library (sctransform)

### Slide 1
brain <- readRDS ("brain_slide1_G1G3_groups.rds")

brain <- SCTransform(brain, vst.flavor = "v2") 

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

brain2 <- SCTransform(brain2, vst.flavor = "v2") 


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



#######
## The viridis coloring

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
ggsave (paste (gene, " plot.pdf", sep=""), pa3, width=8, height=8)



######
## see https://stackoverflow.com/questions/70942728/understanding-color-scales-in-ggplot2
## the original version is based on the viridis palette

p1o <- SpatialFeaturePlot(brain, images=c("slice1"), features = gene)
p2o <- SpatialFeaturePlot(brain, images=c("slice1"), features = gene) + ggplot2::scale_color_viridis_c(option= "rocket")
p1o | p2o

p1 | p1o




######

### The scale_fill_gradient version (more control of the colors, i.e each slide has a common range for the gradient)


add_image <- function (gene) {
max_gene <- max (brain[["SCT"]]$data [row.names (brain[["SCT"]]$data) == gene, ])
midpoint <- max_gene /2


p1 <- SpatialFeaturePlot(brain, images=c("slice1"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p2 <- SpatialFeaturePlot(brain, images=c("slice1.3"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

p3 <- SpatialFeaturePlot(brain, images=c("slice1.2"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p4 <- SpatialFeaturePlot(brain, images=c("slice1.4"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

pa1 <- ggarrange (p1 | p2, nrow=1, labels="G1")
pa2 <- ggarrange (p3 | p4, nrow=1, labels="G3")
pall1 <- ggarrange (pa1, pa2, nrow=2)



rm (p1,p2,p3,p4)

p1 <- SpatialFeaturePlot(brain2, images=c("slice1"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p2 <- SpatialFeaturePlot(brain2, images=c("slice1.3"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

p3 <- SpatialFeaturePlot(brain2, images=c("slice1.2"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p4 <- SpatialFeaturePlot(brain2, images=c("slice1.4"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

pa1 <- ggarrange (p1 | p2, nrow=1, labels="G2")
pa2 <- ggarrange (p3 | p4, nrow=1, labels="G4")
pall2 <- ggarrange (pa1, pa2, nrow=2)


pall3 <- pall1 | pall2
return (pall3)
}

gene <- "Ttr"
pall3 <- add_image (gene)
pall3

ggsave (paste (gene, " version 2 plot.pdf", sep=""), pall3, width=8, height=8)






