library (Seurat)


brain <- readRDS ("brain_G2G1_groups.rds")


########## SCT normalized data
## extract SCT data
counts <- as.matrix (brain[["SCT"]]$data [ ,WhichCells(brain, expression = location == "Hippocampus")])
dim (counts)
# 18827   229


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

boxplot (res$log.fold.change)
abline (h=0)
res.mean <- res


boxplot (apply (counts, 1, mean))
counts <- counts[apply (counts, 1, mean) > 0.2, ]
dim (counts)
# 8568  229

meta <- meta[row.names (meta) %in% colnames (counts), ]
stopifnot (row.names (meta) == colnames (counts))

res.mean <- res.mean[row.names (res.mean) %in% row.names (counts), ]



## For wilcoxon test, see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

counts.s1 <- counts[ ,grepl("G1", meta$group)]
counts.s2 <- counts[ ,grepl("G2", meta$group)]

resl <- list ()
for (i in (1:dim(counts.s1)[1])) {
res <- cbind (data.frame (gene_name= row.names (counts)[i]), data.frame (p.value= wilcox.test(counts.s1[i, ], counts.s2[i ,], exact=FALSE)$p.value))
resl[[i]] <- res
}

res <- do.call ("rbind", resl)
colnames (res)[1] <- "pval.wilcox"
res <- merge (res.mean, res, by="row.names") 
res$padj <- p.adjust (res$pval.wilcox, method="BH")

table (res$padj < 0.05)
#FALSE  TRUE 
# 5624  2944 






