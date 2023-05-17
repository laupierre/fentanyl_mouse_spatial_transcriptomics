## This was taken from: deg_lognorm_brain_normalization_hippocampus_G2vsG1.R

library (Seurat)
library (scuttle)
library(scran)
library (scater)
library (SingleCellExperiment)
library (openxlsx)

brain <- readRDS ("brain_slide2_G2G4_groups.rds")
brain

sce <- as.SingleCellExperiment(brain)

# add a low expression gene filtering
ave.counts <- rowMeans(counts(sce)) 

# look at a chosen log10 threshold, here 0.05, with log10 (0.05) = -1.3
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab= expression(Log[10]~"average count"))
abline(v=log10(0.05), col="blue", lwd=2, lty=2)

sce <- sce[ave.counts >= 0.05, ]
sce


clusters <- quickCluster(sce)

sce <- computePooledFactors(sce, clusters=clusters)
summary(sizeFactors(sce))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5281  0.8261  0.9878  1.0000  1.1882  1.4870 

# normalization of the entire brain
sce <- logNormCounts(sce)
counts <- logcounts (sce)

meta <- colData (sce)
idx <- match (row.names (meta), colnames (counts))
meta <- meta[idx, ]
stopifnot (row.names (meta) == colnames (counts))

table (meta$group)

## See the centering of fold changes in all the cells
counts.s <- counts[ ,grepl("G2", meta$group)]
mean.s1 <- apply (counts.s, 1, mean)
counts.s <- counts[ ,grepl("G4", meta$group)]
mean.s2 <- apply (counts.s, 1, mean)

res <- cbind (data.frame (mean.s2), data.frame (mean.s1))
res <-  cbind (res, data.frame (log.fold.change= res$mean.s2 - res$mean.s1))
colnames (res)[3] <- "log.fold.change"

boxplot (res$log.fold.change)
abline (h=0)
res.mean <- res



## Extract the cells of interest, i.e hippocampus

meta.s <- meta[meta$location == "Hippocampus", ]

counts <- counts[ ,colnames (counts) %in% row.names (meta.s)]
idx <- match (row.names (meta.s), colnames (counts))
meta.s <- meta.s[idx, ]
stopifnot (row.names (meta.s) == colnames (counts))


## See the centering of fold changes in the selected cells only 
counts.s <- counts[ ,grepl("G2", meta.s$group)]
mean.s1 <- apply (counts.s, 1, mean)
counts.s <- counts[ ,grepl("G4", meta.s$group)]
mean.s2 <- apply (counts.s, 1, mean)

res <- cbind (data.frame (mean.s2), data.frame (mean.s1))
res <-  cbind (res, data.frame (log.fold.change= res$mean.s2 - res$mean.s1))
colnames (res)[3] <- "log.fold.change"

boxplot (res$log.fold.change)
abline (h=0)
res.mean <- res


## For wilcoxon test, see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

counts.s1 <- counts[ ,grepl("G4", meta.s$group)]
counts.s2 <- counts[ ,grepl("G2", meta.s$group)]

resl <- list ()
for (i in (1:dim(counts.s1)[1])) {
if (i %% 500 == 0) {
print (i)
}
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
#12158    36


## Gene annotation

library('org.Mm.eg.db')

#columns(org.Mm.eg.db)
symbols <- res$gene_name
res1a <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (res$gene_name, names (res1a))
res$Description <- as.vector (res1a) [idx]
res <- res[order (res$padj), ]

write.xlsx (res, "table 6. hippocampus_G4vsG2_selected_cells_brain_normalization_wilcoxon_analysis.xlsx", rowNames=F)





