library (openxlsx)

## bulk RNA-Seq
bulk <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")
colnames (bulk)[2] <- "gene_name"

## Pseudobulk RNA-Seq
pseudo <- read.xlsx ("pseudobulk_hippocampus_selected_cells.xlsx")
pseudo <- pseudo[ ,c("gene_name", "logFC", "adj.P.Val", "AveExpr")]

## Single-cell (from spatial transcriptomics)
sc <- read.xlsx ("hippocampus_selected_cells_wilcoxon_analysis.xlsx")
# invert the log fold changes
sc$avg_logFC <- -1*sc$avg_logFC
sc <- sc[ ,c("gene_name", "avg_logFC", "p_val_adj", "mean", "Description")]

par (mfrow= c(3,1))
comp1 <- merge (sc, pseudo, by="gene_name")
plot (comp1$avg_logFC, comp1$logFC, xlab="log fold changes wilcoxon (spatial)", ylab="log fold changes limma (pseudobulk)", main="Comparison spatial vs pseudobulk transcriptomics")
abline (0,1, col="red")
abline (h=0)
abline (v=0)

comp2 <- merge (sc, bulk, by="gene_name")
plot (comp2$avg_logFC, comp2$Log2.Fold.Change, xlab="log fold changes wilcoxon (spatial)", ylab="log fold changes Chicago (bulk)", main="Comparison spatial vs bulk transcriptomics")
abline (0,1, col="red")
abline (h=0)
abline (v=0)

comp3 <- merge (pseudo, bulk, by="gene_name")
plot (comp2$logFC, comp2$Log2.Fold.Change, xlab="log fold changes limma (pseudobulk)", ylab="log fold changes Chicago (bulk)", main="Comparison pseudobulk vs bulk transcriptomics")
abline (0,1, col="red")
abline (h=0)
abline (v=0)


