### SCT is not working well because significant genes are not symmetrical 

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


#boxplot (apply (counts, 1, mean))
## The expression normalization is post SCT test (not usual)!! 
#counts <- counts[apply (counts, 1, mean) > 0.2, ]
#dim (counts)
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



#########
## Sanity check. comparison with bulk RNA-Seq and normalization on selected cells

library (openxlsx)

## bulk RNA-Seq
bulk <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")
colnames (bulk)[2] <- "gene_name"

## single cell (normalized on selected hippocampal cells)
sc <- read.xlsx ("hippocampus_selected_cells_wilcoxon_analysis.xlsx")
# invert the log fold changes
sc$avg_logFC <- -1*sc$avg_logFC
sc <- sc[ ,c("gene_name", "avg_logFC", "p_val_adj", "mean", "Description")]


comp1 <- merge (res, bulk, by="gene_name") 
plot (comp1$log.fold.change, comp1$Log2.Fold.Change, xlab="log fold changes wilcoxon (spatial)", ylab="log fold changes wald (bulk)", main="Comparison spatial vs bulk transcriptomics",
      xlim=c(-3,3), ylim=c(-3,3))
abline (0,1, col="red")
abline (h=0)
abline (v=0)



### SCT log fold changes are compressed relative to hippocampal normalization

comp2 <- merge (res, sc, by="gene_name") 
plot (comp2$log.fold.change, comp2$avg_logFC, xlab="log fold changes wilcoxon (brain norm)", ylab="log fold changes wilcoxon (hippocampus norm)", main="Comparison spatial vs bulk transcriptomics",
      xlim=c(-3,3), ylim=c(-3,3), col=ifelse (comp2$padj < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)








