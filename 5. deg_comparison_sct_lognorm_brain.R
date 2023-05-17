######### Comparison SCT and lognorm brain normalization

library (openxlsx)

sct <- read.xlsx ("table 7. hippocampus_G2vsG1_selected_cells_sct_brain_normalization_wilcoxon_analysis.xlsx")
log.n <- read.xlsx ("table 4. hippocampus_G2vsG1_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (sct, log.n, by="gene_name") 

pdf ("figure 1. Comparison SCT and log normalization methods G2 vs G1.pdf")
plot (comp2$log.fold.change.x, comp2$log.fold.change.y, xlab="logFC wilcoxon (SCT normalization)", ylab="logFC (lognorm norm)", main="Comparison SCT and log normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj.x < 0.05 | comp2$padj.y < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()


sct <- read.xlsx ("table 8. hippocampus_G3vsG1_selected_cells_sct_brain_normalization_wilcoxon_analysis.xlsx")
log.n <- read.xlsx ("table 5. hippocampus_G3vsG1_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (sct, log.n, by="gene_name") 

pdf ("figure 2. Comparison SCT and log normalization methods G3 vs G1.pdf")
plot (comp2$log.fold.change.x, comp2$log.fold.change.y, xlab="logFC wilcoxon (SCT normalization)", ylab="logFC (lognorm norm)", main="Comparison SCT and log normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj.x < 0.05 | comp2$padj.y < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()


sct <- read.xlsx ("table 9. hippocampus_G4vsG2_selected_cells_sct_brain_normalization_wilcoxon_analysis.xlsx")
log.n <- read.xlsx ("table 6. hippocampus_G4vsG2_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (sct, log.n, by="gene_name") 

pdf ("figure 3. Comparison SCT and log normalization methods G4 vs G2.pdf")
plot (comp2$log.fold.change.x, comp2$log.fold.change.y, xlab="logFC wilcoxon (SCT normalization)", ylab="logFC (lognorm norm)", main="Comparison SCT and log normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj.x < 0.05 | comp2$padj.y < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()



######### Comparison lognorm on hippocampal cells and lognorm brain normalization
rm (list= ls ())

log.n <- read.xlsx ("table 4. hippocampus_G2vsG1_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")
## single cell (normalized on selected hippocampal cells and not on the entire brain)
log.n.sc <- read.xlsx ("table 1. hippocampus_G2vsG1_selected_cells_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (log.n, log.n.sc, by="gene_name") 

pdf ("figure 4. Comparison lognorm entire brain and selected cells methods G2 vs G1.pdf")
plot (comp2$log.fold.change, comp2$avg_logFC, xlab="logFC wilcoxon (lognorm entire brain)", ylab="logFC (lognorm on selected cells)", main="Comparison entire brain and selected cells normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj < 0.05 | comp2$p_val_adj < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()


log.n <- read.xlsx ("table 5. hippocampus_G3vsG1_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")
## single cell (normalized on selected hippocampal cells and not on the entire brain)
log.n.sc <- read.xlsx ("table 2. hippocampus_G3vsG1_selected_cells_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (log.n, log.n.sc, by="gene_name") 

pdf ("figure 5. Comparison lognorm entire brain and selected cells methods G3 vs G1.pdf")
plot (comp2$log.fold.change, comp2$avg_logFC, xlab="logFC wilcoxon (lognorm entire brain)", ylab="logFC (lognorm on selected cells)", main="Comparison entire brain and selected cells normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj < 0.05 | comp2$p_val_adj < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()



log.n <- read.xlsx ("table 6. hippocampus_G4vsG2_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")
## single cell (normalized on selected hippocampal cells and not on the entire brain)
log.n.sc <- read.xlsx ("table 3. hippocampus_G4vsG2_selected_cells_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (log.n, log.n.sc, by="gene_name") 

pdf ("figure 6. Comparison lognorm entire brain and selected cells methods G4 vs G2.pdf")
plot (comp2$log.fold.change, comp2$avg_logFC, xlab="logFC wilcoxon (lognorm entire brain)", ylab="logFC (lognorm on selected cells)", main="Comparison entire brain and selected cells normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj < 0.05 | comp2$p_val_adj < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()











