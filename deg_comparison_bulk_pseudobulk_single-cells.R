library (opnxslx)

## bulk RNA-Seq
bulk <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")

## Pseudobulk RNA-Seq
pseudo <- read.xlsx ("pseudobulk_hippocampus_selected_cells.xlsx")

## Single-cell (from spatial transcriptomics)
sc <- read.xlsx ("hippocampus_selected_cells_wilcoxon_analysis.xlsx")

