## See the scuttle library: https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/norm.html
## This method keeps the symmetry in comparison to the SCT normalization (see step 4)

library (Seurat)
library (scuttle)
library(scran)
library (scater)
library (SingleCellExperiment)
library (openxlsx)

brain <- readRDS ("brain_G2G1_groups.rds")
brain

#counts <- brain[["Spatial"]]$counts
#colData <- brain@meta.data
#sce <- SingleCellExperiment(assays = list(counts = counts, colData=colData))

# The counts slot of the "SCT" assay is supposed to be adjusted/corrected counts derived by reverse-transforming the Pearson residuals
# The counts slot of the "RNA" assay, on the other hand, would just have the raw counts read from the UMI data
# So, they can be different, and in this case, the SCT-normalized corrected counts for some of the genes end up increasing post normalization. 
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
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4482  0.8359  0.9831  1.0000  1.1656  1.6327 


# normalization of the entire brain
sce <- logNormCounts(sce)
counts <- logcounts (sce)

meta <- colData (sce)
idx <- match (row.names (meta), colnames (counts))
meta <- meta[idx, ]
stopifnot (row.names (meta) == colnames (counts))

## See the centering of fold changes in all the cells
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


## Extract the cells of interest, i.e hippocampus

meta.s <- meta[meta$location == "Hippocampus", ]

counts <- counts[ ,colnames (counts) %in% row.names (meta.s)]
idx <- match (row.names (meta.s), colnames (counts))
meta.s <- meta.s[idx, ]
stopifnot (row.names (meta.s) == colnames (counts))


## See the centering of fold changes in the selected cells only 
counts.s <- counts[ ,grepl("G1", meta.s$group)]
mean.s1 <- apply (counts.s, 1, mean)
counts.s <- counts[ ,grepl("G2", meta.s$group)]
mean.s2 <- apply (counts.s, 1, mean)

res <- cbind (data.frame (mean.s2), data.frame (mean.s1))
res <-  cbind (res, data.frame (log.fold.change= res$mean.s2 - res$mean.s1))
colnames (res)[3] <- "log.fold.change"

boxplot (res$log.fold.change)
abline (h=0)
res.mean <- res



## For wilcoxon test, see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

counts.s1 <- counts[ ,grepl("G1", meta.s$group)]
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
#12465    91  

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

write.xlsx (res, "table 4. hippocampus_selected_cells_brain_normalization_wilcoxon_analysis.xlsx", rowNames=F)



## Sanity check (with bulk RNA-Seq)

rnaseq <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")
rnaseq <- rnaseq[ ,c("Gene.Symbol", "Log2.Fold.Change", "FDR.Adj.p.Value", "Mean")]

res2 <- merge (res, rnaseq, by.x="gene_name", by.y="Gene.Symbol")

# write.xlsx (res2, "hippocampus_selected_cells_comparison_rnaseq_vs_wilcoxon_analysis_brain_normalization.xlsx", rowNames=F)

plot (res2$log.fold.change, res2$Log2.Fold.Change, xlab="Wilcoxon brain normalization", ylab="Chicago RNA-Seq", col=ifelse (res2$padj < 0.05, "blue", "black"),
      xlim=c(-3,3), ylim=c(-3,3))
abline (h=0)
abline (v=0)



## single cell (normalized on selected hippocampal cells and not on the entire brain)
sc <- read.xlsx ("hippocampus_selected_cells_wilcoxon_analysis.xlsx")
# invert the log fold changes
sc$avg_logFC <- -1*sc$avg_logFC
sc <- sc[ ,c("gene_name", "avg_logFC", "p_val_adj", "mean", "Description")]

comp1 <- merge (res, sc, by="gene_name") 
plot (comp1$log.fold.change, comp1$avg_logFC, xlab="wilcoxon brain normalization", ylab="wilcoxon hippocampus normalization", main="Comparison spatial vs bulk transcriptomics",
      col=ifelse (comp1$padj < 0.05, "blue", "black"),
      xlim=c(-3,3), ylim=c(-3,3))
abline (0,1, col="red")
abline (h=0)
abline (v=0)











