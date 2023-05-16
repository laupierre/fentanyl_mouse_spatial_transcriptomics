## See the scuttle library: https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/norm.html

library (Seurat)
library (scuttle)
library(scran)
library (scater)
library (SingleCellExperiment)


brain <- readRDS ("brain_G2G1_groups.rds")
brain

#counts <- brain[["Spatial"]]$counts
#colData <- brain@meta.data
#sce <- SingleCellExperiment(assays = list(counts = counts, colData=colData))

# The counts slot of the "SCT" assay is supposed to be adjusted/corrected counts derived by reverse-transforming the Pearson residuals
# The counts slot of the "RNA" assay, on the other hand, would just have the raw counts read from the UMI data
# So, they can be different, and I think in your case, the SCT-normalized corrected counts for some of the genes you mentioned end up increasing post normalization. 
# Essentially, the SCT@counts may be viewed as the normalized versions of the RNA@counts based on the SCT technique

#counts <- brain[["Spatial"]]$counts
#head (counts [ ,1:10])

#counts <- brain[["SCT"]]$counts
#head (counts [ ,1:10])

# this parses the sequencing depth corrected counts from SCT, ie $counts
sce <- as.SingleCellExperiment(brain)
head (counts(sce) [ ,1:10])

# add a low expression gene filtering
ave.counts <- rowMeans(counts(sce)) 

# look at a chosen log10 threshold, here 0.05, with log10 (0.05) = -1.3
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab= expression(Log[10]~"average count"))
abline(v=log10(0.05), col="blue", lwd=2, lty=2)

sce <- sce[ave.counts >= 0.05, ]
sce


# or select genes that have non-zero counts in at least n cells (ie 1000)

#numcells <- nexprs(sce, byrow=TRUE)
#alt.keep <- numcells >= 1000
#sum(alt.keep)
#smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), ylab="Number of expressing cells")


clusters <- quickCluster(sce)

sce <- computePooledFactors(sce, clusters=clusters)
summary(sizeFactors(sce))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4400  0.8341  0.9848  1.0000  1.1655  1.6719


# normalization of the entire brain
sce <- logNormCounts(sce)
counts <- logcounts (sce)

meta <- colData (sce)
#meta.s <- meta [meta$location == "Hippocampus", ]
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
#14007  2637 

boxplot (res$log.fold.change)
abline (h=0)












