library (opnxslx)

## bulk RNA-Seq
bulk <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")

## Pseudobulk RNA-Seq
pseudo <- read.xlsx ("pseudobulk_hippocampus_selected_cells.xlsx")
pseudo <- pseudo[ ,c("gene_name", "logFC", "adj.P.Val", "AveExpr")]

## Single-cell (from spatial transcriptomics)
sc <- read.xlsx ("hippocampus_selected_cells_wilcoxon_analysis.xlsx")
# invert the log fold changes
sc$avg_logFC <- -1*sc$avg_logFC
sc <- sc[ ,c("gene_name", "avg_logFC", "p_val_adj", "mean", "Description")]

comp1 <- merge (sc, pseudo, by="gene_name")
plot (comp1$avg_logFC, comp1$logFC, xlab="wilcoxon", "limma")
abline (h=0)
abline (v=0)

