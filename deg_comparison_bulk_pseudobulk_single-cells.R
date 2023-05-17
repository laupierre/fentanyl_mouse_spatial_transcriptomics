library (openxlsx)

## bulk RNA-Seq
bulk <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")
colnames (bulk)[2] <- "gene_name"

## Pseudobulk
pseudo <- read.xlsx ("pseudobulk_hippocampus_lognorm_selected_cells.xlsx")
pseudo <- pseudo[ ,c("gene_name", "logFC", "adj.P.Val", "AveExpr")]

## Single-cell (from lognorm spatial transcriptomics)
sc <- read.xlsx ("table 1. hippocampus_G2vsG1_selected_cells_normalization_wilcoxon_analysis.xlsx")
# invert the log fold changes
#sc$avg_logFC <- -1*sc$avg_logFC
sc <- sc[ ,c("gene_name", "avg_logFC", "p_val_adj", "mean", "Description")]


pdf ("Comparison between bulk, pseudobulk and spatial RNA-Seq version 2.pdf")
par (mfrow= c(3,1))
comp1 <- merge (sc, pseudo, by="gene_name")
plot (comp1$avg_logFC, comp1$logFC, xlab="logFC wilcoxon (lognorm selected cells)", ylab="logFC limma (pseudobulk)", main="Comparison spatial vs pseudobulk transcriptomics",
      xlim=c(-3,3), ylim=c(-3,3))
abline (0,1, col="red")
abline (h=0)
abline (v=0)

comp2 <- merge (sc, bulk, by="gene_name")
plot (comp2$avg_logFC, comp2$Log2.Fold.Change, xlab="logFC wilcoxon (lognorm selected cells)", ylab="logFC Bulk RNA-Seq", main="Comparison spatial vs bulk transcriptomics",
     xlim=c(-3,3), ylim=c(-3,3))
abline (0,1, col="red")
abline (h=0)
abline (v=0)

comp3 <- merge (pseudo, bulk, by="gene_name")
plot (comp3$logFC, comp3$Log2.Fold.Change, xlab="logFC limma (pseudobulk)", ylab="logfold Bulk RNA-Seq", main="Comparison pseudobulk vs bulk transcriptomics",
     xlim=c(-3,3), ylim=c(-3,3))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()


# Internal table comparison

comp1 <- merge (sc, pseudo, by="gene_name")
comp2 <- merge (comp1, bulk, by="gene_name")
comp2$trend <- ifelse (comp2$avg_logFC >0 & comp2$logFC >0 & comp2$Log2.Fold.Change >0, "more_G2", ifelse (comp2$avg_logFC <0 & comp2$logFC <0 & comp2$Log2.Fold.Change <0, "less_G2", "No"))
table (comp2$trend)
# less_G2 more_G2      No 
# 1983    1745    4209

comp2 <- comp2[order (comp2$avg_logFC), ]
colnames (comp2)[1:5] <- paste (colnames (comp2)[1:5], "wilcoxon", sep="-")
colnames (comp2)[6:8] <- paste (colnames (comp2)[6:8], "pseudobulk", sep="-")
colnames (comp2)[9:34] <- paste (colnames (comp2)[9:34], "wald_bulk", sep="-")
head (comp2)

write.xlsx (comp2, "Comparison between bulk, pseudobulk and spatial RNA-Seq.xlsx", rowNames=F)







