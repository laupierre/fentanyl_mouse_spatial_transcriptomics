library (Seurat)
library (ggpubr)
library (sctransform)
library (RColorBrewer)

### Slide 1
brain <- readRDS ("brain_slide1_G1G3_groups.rds")
### Slide 2
brain2 <- readRDS ("brain_slide2_G2G4_groups.rds")

colors <- c(Hippocampus = brewer.pal(3, "Dark2")[2], unknown = "grey")

p1 <- SpatialDimPlot(brain, cols= colors, group.by="location", images = "slice1") + ggplot2::theme(legend.position = "none")
p2 <- SpatialDimPlot(brain, cols= colors, group.by="location", images = "slice1.3") + ggplot2::theme(legend.position = "none")
p3 <- SpatialDimPlot(brain, cols= colors, group.by="location", images = "slice1.2") + ggplot2::theme(legend.position = "none")
p4 <- SpatialDimPlot(brain, cols= colors, group.by="location", images = "slice1.4") + ggplot2::theme(legend.position = "none")

pa1 <- ggarrange (p1 | p2, nrow=1, labels="G1")
pa2 <- ggarrange (p3 | p4, nrow=1, labels="G3")
pall1 <- ggarrange (pa1, pa2, nrow=2)


rm (p1,p2,p3,p4)

p1 <- SpatialDimPlot(brain2, cols= colors, group.by="location", images = "slice1") + ggplot2::theme(legend.position = "none")
p2 <- SpatialDimPlot(brain2, cols= colors, group.by="location", images = "slice1.3") + ggplot2::theme(legend.position = "none")
p3 <- SpatialDimPlot(brain2, cols= colors, group.by="location", images = "slice1.2") + ggplot2::theme(legend.position = "none")
p4 <- SpatialDimPlot(brain2, cols= colors, group.by="location", images = "slice1.4") + ggplot2::theme(legend.position = "none")

pa1 <- ggarrange (p1 | p2, nrow=1, labels="G2")
pa2 <- ggarrange (p3 | p4, nrow=1, labels="G4")
pall2 <- ggarrange (pa1, pa2, nrow=2)

ggsave ("selected cells highlight 1 plot.pdf", sep=""), pall1, width=8, height=8)
ggsave ("selected cells highlight 2 plot.pdf", sep=""), pall2, width=8, height=8)


