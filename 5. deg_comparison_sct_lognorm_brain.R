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

