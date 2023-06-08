## This is the screening of the Hippo-Seq data
## To increase the speed of screening, we output two slices only (out of 4)

library (Seurat)
library (ggpubr)
library (sctransform)


### Pseudo-slide

brain <- readRDS ("brain_G2G1_groups.rds")

brain <- SCTransform(brain, vst.flavor = "v2") 

Images (brain)
# "slice1"   "slice1.2" "slice1.3" "slice1.4"
# G2-2C, G2-2A, G1-1C, G1-1A

add_image <- function (gene) {
max_gene <- max (brain[["SCT"]]$data [row.names (brain[["SCT"]]$data) == gene, ])
midpoint <- max_gene /2

p1 <- SpatialFeaturePlot(brain, images=c("slice1"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))
p2 <- SpatialFeaturePlot(brain, images=c("slice1.3"), features = gene) + ggplot2::scale_fill_gradient2(midpoint = midpoint, low="blue", mid="white", high="red", limits = c(0,max_gene),
						 breaks = round (seq(0, max_gene, length.out = 6), digits=1))

p3 <- ggarrange (p1 | p2, nrow=1)
return (p3)
}

## Hippo-Seq data
genes <- readLines ("hipposeq_genes.txt")

for (i in (1:length (genes))) {
print (i)
gene <- genes[i]
tryCatch (
{
p <- add_image (gene)
print (p)
invisible(readline(prompt="Press [enter] to continue"))
}, error=function(e) {cat("error\n")}, warning=function (w){cat("warning\n")}
)
}


