### SCT is not working well because there are many significant genes that are not symmetrical 
### See info: https://github.com/satijalab/seurat/issues/5009#issuecomment-907426252
### If your library sizes are comparable across the different slides, you can compare the data slots of SCT assay. 
### If the median sequencing depths are similar, the data slot of SCT assay is comparable.
### If you do not expect huge slide to slide variation in terms of technical noise, you would get very similar results if you run SCTransform after merging. 
### The idea behind splitting and then running SCTransform is to enable it to learn a dataset-specific model of technical noise (which could be very similar across samples in most cases)

### See the SCTtransform version 2, https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

# install.packages("sctransform", dependencies=TRUE)

library (sctransform)
library (Seurat)
library (openxlsx)

brain <- readRDS ("brain_G2G1_groups.rds")

brain <- SCTransform(brain, vst.flavor = "v2") 
 
 
####### SCT normalized data
## extract the SCT hippocampus normalized data directly
counts <- as.matrix (brain[["SCT"]]$data [ ,WhichCells(brain, expression = location == "Hippocampus")])
dim (counts)
# 18827   229

## There's no low expression filtering here!
## manual calculation of SCT log fold changes and wilcoxon test of SCT values

meta <- brain [ ,WhichCells(brain, expression = location == "Hippocampus")]@meta.data
idx <- match (row.names (meta), colnames (counts))
meta <- meta[idx, ]
stopifnot (row.names (meta) == colnames (counts))

counts.s <- counts[ ,grepl("G1", meta$group)]
mean.s1 <- apply (counts.s, 1, mean)
counts.s <- counts[ ,grepl("G2", meta$group)]
mean.s2 <- apply (counts.s, 1, mean)

res <- cbind (data.frame (mean.s2), data.frame (mean.s1))
res <-  cbind (res, data.frame (log.fold.change= res$mean.s2 - res$mean.s1))
colnames (res)[3] <- "log.fold.change"

res <- res[res$mean.s2 != 0 & res$mean.s1 != 0, ]
#table (res$mean.s2 == 0 & res$mean.s1 == 0)

boxplot (res$log.fold.change)
abline (h=0)
res.mean <- res



## For wilcoxon test, see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

counts.s1 <- counts[ ,grepl("G1", meta$group)]
counts.s2 <- counts[ ,grepl("G2", meta$group)]

resl <- list ()
for (i in (1:dim(counts.s1)[1])) {
res <- cbind (data.frame (gene_name= row.names (counts)[i]), data.frame (p.value= wilcox.test(counts.s1[i, ], counts.s2[i ,], exact=FALSE)$p.value))
resl[[i]] <- res
}

res <- do.call ("rbind", resl)
colnames (res)[2] <- "pval.wilcox"
res <- merge (res.mean, res, by.x="row.names", by.y="gene_name") 
res$padj <- p.adjust (res$pval.wilcox, method="BH")
res <- res[order (res$padj), ]
colnames (res)[1] <- "gene_name"
head (res)

table (res$padj < 0.05)
#FALSE  TRUE 
#12576   452   

boxplot (res$log.fold.change)
abline (h=0)


## Gene annotation

library('org.Mm.eg.db')

#columns(org.Mm.eg.db)
symbols <- res$gene_name
res1a <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (res$gene_name, names (res1a))
res$Description <- as.vector (res1a) [idx]
res <- res[order (res$padj), ]

write.xlsx (res, "table 7. hippocampus_G2vsG1_selected_cells_sct_brain_normalization_wilcoxon_analysis.xlsx", rowNames=F)












