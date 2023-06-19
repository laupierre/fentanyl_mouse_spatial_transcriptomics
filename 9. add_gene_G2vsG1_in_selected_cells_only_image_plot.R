library (Seurat)
library (ggpubr)
library (sctransform)


### Pseudo-slide

brain <- readRDS ("brain_G2G1_groups.rds")

brain <- SCTransform(brain, vst.flavor = "v2") 

Images (brain)
# "slice1"   "slice1.2" "slice1.3" "slice1.4"
# G2-2C, G2-2A, G1-1C, G1-1A

# subset to the selected hippocampal spots only

Idents (brain) <- "location"
brain <- subset (brain, idents = c("Hippocampus"))



## without cropping

add_image <- function (gene) {
max_gene <- max (brain[["SCT"]]$data [row.names (brain[["SCT"]]$data) == gene, ])
midpoint <- max_gene /2


p1 <- SpatialFeaturePlot(brain, images=c("slice1"), crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p2 <- SpatialFeaturePlot(brain, images=c("slice1.2"), crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

p3 <- SpatialFeaturePlot(brain, images=c("slice1.3"), crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p4 <- SpatialFeaturePlot(brain, images=c("slice1.4"), crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

pa1 <- ggarrange (p2 | p1, nrow=1, labels="G2")
pa2 <- ggarrange (p4 | p3, nrow=1, labels="G1")
pall1 <- ggarrange (pa1, pa2, nrow=2)

return (pall1)
}

gene <- "Plppr4"
pall3 <- add_image (gene)
pall3

ggsave (paste (gene, " version 3 plot.pdf", sep=""), pall3, width=8, height=8)



## with cropping, the figure is stitched, see https://github.com/satijalab/seurat/issues/5212

add_image <- function (gene) {
max_gene <- max (brain[["SCT"]]$data [row.names (brain[["SCT"]]$data) == gene, ])
midpoint <- max_gene /2


p1 <- SpatialFeaturePlot(brain, images=c("slice1"), crop=TRUE, pt.size.factor = 4, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p2 <- SpatialFeaturePlot(brain, images=c("slice1.2"), crop=TRUE, pt.size.factor = 4, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

p3 <- SpatialFeaturePlot(brain, images=c("slice1.3"), crop=TRUE, pt.size.factor = 4, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p4 <- SpatialFeaturePlot(brain, images=c("slice1.4"), crop=TRUE, pt.size.factor = 4, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

pa1 <- ggarrange (p2 | p1, nrow=1, labels="G2")
pa2 <- ggarrange (p4 | p3, nrow=1, labels="G1")
pall1 <- ggarrange (pa1, pa2, nrow=2)

return (pall1)
}

gene <- "Plppr4"
pall3 <- add_image (gene)
pall3

ggsave (paste (gene, " version 4 plot.pdf", sep=""), pall3, width=8, height=8)




## without cropping and rotating the image to be on correct orientation

add_image <- function (gene) {
max_gene <- max (brain[["SCT"]]$data [row.names (brain[["SCT"]]$data) == gene, ])
midpoint <- max_gene /2

p1 <- SpatialFeaturePlot(brain, images=c("slice1"), pt.size.factor = 1, crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p2 <- SpatialFeaturePlot(brain, images=c("slice1.2"), pt.size.factor = 1, crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p3 <- SpatialFeaturePlot(brain, images=c("slice1.3"), pt.size.factor = 1, crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p4 <- SpatialFeaturePlot(brain, images=c("slice1.4"), pt.size.factor = 1, crop=FALSE, alpha = c(0.8, 1), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene), 
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

pa1 <- ggarrange (p2 | p1, nrow=1, labels="")
pa2 <- ggarrange (p4 | p3, nrow=1, labels="")
pall1 <- ggarrange (pa1, pa2, nrow=2)

return (pall1)
}


gene <- "Ttr"
pall3 <- add_image (gene)


## Rotate the grid

library (grid)
library (gridExtra)
library (ggpubr)

grid.newpage() 
pushViewport(viewport(angle=-90, width = unit(8, "inches"), height = unit(8, "inches")))  
grid.draw(ggplot_gtable(ggplot_build(pall3)))


pdf (paste (gene, "version 5 rotated plot.pdf"), height = 8, width = 8)
grid.newpage() 
pushViewport(viewport(angle=-90, width = unit(8, "inches"), height = unit(8, "inches")))  
grid.draw(ggplot_gtable(ggplot_build(pall3)))
dev.off ()



# subset the selected hippocampal spots only

Idents (brain) <- "location"
brain <- subset (brain, idents = c("Hippocampus"))

gene <- "Ttr"
pall3 <- add_image (gene)


## Rotate the grid

library (grid)
library (gridExtra)
library (ggpubr)

grid.newpage() 
pushViewport(viewport(angle=-90, width = unit(8, "inches"), height = unit(8, "inches")))  
grid.draw(ggplot_gtable(ggplot_build(pall3)))


pdf (paste (gene, "version 6 rotated plot.pdf"), height = 8, width = 8)
grid.newpage() 
pushViewport(viewport(angle=-90, width = unit(8, "inches"), height = unit(8, "inches")))  
grid.draw(ggplot_gtable(ggplot_build(pall3)))
dev.off ()






