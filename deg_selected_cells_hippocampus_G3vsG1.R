library (Libra)
library (Seurat)
library (openxlsx)

brain <- readRDS ("brain_slide1_G1G3_groups.rds")

# raw counts
counts <- as.matrix (brain[["SCT"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])
dim (counts)
# 18966   252

counts <- counts[apply (counts, 1, mean) > 0.05, ]
dim (counts)
# 12691   252


# for backward compatibility
counts <- as (counts, "sparseMatrix")
seurat.ss2 <- CreateSeuratObject(counts)

meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == "Hippocampus"), ]
# we need cell_type (eg hippocampus), replicate (eg mouse), and label (eg treatment)
meta.ss2$replicate <- gsub (".*-", "", meta.ss2$group)
meta.ss2$label <- gsub ("-.*", "", meta.ss2$group)
meta.ss2$cell_type <- "Hippocampus"

meta.ss2 <- meta.ss2[row.names (meta.ss2) %in% colnames (counts), ]
idx <- match (row.names (meta.ss2), colnames (counts))
meta.ss2 <- meta.ss2[idx, ]
stopifnot (row.names (meta.ss2) == colnames (counts))

seurat.ss2@meta.data <- meta.ss2
table (seurat.ss2@meta.data$label)   # It will compare G1 vs G3
# G1  G3 
#120 132 

mymean <- data.frame (mean= apply (counts, 1, mean))


res <- run_de(seurat.ss2, de_method = 'wilcox', de_family= "singlecell")
res <- data.frame (res)
row.names (res) <- res$gene
res <- merge (res, mymean, by="row.names")
res <- res[order (res$p_val_adj), ]
row.names (res) <- res$gene
colnames (res)[1] <- "gene_name"
res$avg_logFC <- -1 * res$avg_logFC


## Annotation

library('org.Mm.eg.db')

#columns(org.Mm.eg.db)
symbols <- row.names (res)
res1a <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (row.names (res), names (res1a))
res$Description <- as.vector (res1a) [idx]
res <- res[order (res$p_val_adj), ]
res <- res[ ,-which (colnames (res) == "de_family")]

write.xlsx (res, "hippocampus_G3vsG1_selected_cells_wilcoxon_analysis.xlsx", rowNames=F)



## Sanity check (comparison with brain normalization)

brain.n <- read.xlsx ("hippocampus_G3vsG1_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (res, brain.n, by="gene_name") 

pdf ("Comparison hippocampus and normalization methods G3 vs G1.pdf")
plot (comp2$log.fold.change, comp2$avg_logFC, xlab="log fold changes wilcoxon (brain norm)", ylab="log fold changes wilcoxon (hippocampus norm)", main="Comparison hippocampus and normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()



