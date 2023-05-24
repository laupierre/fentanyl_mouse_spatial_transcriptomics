library (Libra)
library (Seurat)
library (openxlsx)

brain <- readRDS ("brain_slide2_G2G4_groups.rds")

# raw counts
counts <- as.matrix (brain[["SCT"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])
dim (counts)
# 18637   211

counts <- counts[apply (counts, 1, mean) > 0.05, ]
dim (counts)
# 11907   211



#### Optional: Dropout removal

meta <- brain@meta.data
meta.g1 <- meta[grep ("G2", meta$group), ]
meta.g1 <- meta.g1 [row.names (meta.g1) %in% colnames (counts), ]

meta.g2 <- meta[grep ("G4", meta$group), ]
meta.g2 <- meta.g2 [row.names (meta.g2) %in% colnames (counts), ]

counts.g1 <- counts[ ,colnames (counts) %in% row.names (meta.g1)]
counts.g2 <- counts[ ,colnames (counts) %in% row.names (meta.g2)]


# Percentage of cells with zeros
prop1 <- data.frame (prop1= apply (counts.g1, 1, function (x) {sum (x == 0)}))
prop2 <- data.frame (prop2= apply (counts.g2, 1, function (x) {sum (x == 0)}))
prop <- cbind (prop1/dim (counts.g1)[2] *100, prop2/dim (counts.g2)[2] *100)
head (prop)


## see: https://www.jdatalab.com/data_science_and_data_mining/2017/01/30/data-binning-plot.html

library (ggplot2)
library (ggpubr)

p1 <- ggplot(data = prop, mapping = aes(x=prop1)) + 
  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7, bins=20)  +
  geom_density() + xlab ("Percentage of zeros") + ylab ("Density of expressed genes") +  geom_vline(xintercept=75, linetype="dashed", color = "red") +
  ggtitle ("Percentage of zeros in the G2 group")
  
p2 <- ggplot(data = prop, mapping = aes(x=prop2)) + 
  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7, bins=20)  +
  geom_density() + xlab ("Percentage of zeros") + ylab ("Density of expressed genes") + geom_vline(xintercept=75, linetype="dashed", color = "red") +
  ggtitle ("Percentage of zeros in the G4 group") 
  
p3 <- ggarrange (p1, p2, nrow=1)
p3
#ggsave ("figure 12. Percentage of dropouts in G2G4 groups.pdf", p3, width=8, height=8)

# Proportion of cells expressing a gene
prop <- 100- prop

# select genes with at least 25% of cells expressing the gene in any group
prop$drop <- apply (prop, 1, function (x) {any (x > 25)})
prop <- prop[prop$drop == TRUE, ]

counts <- counts[row.names (counts) %in% row.names (prop), ]
dim (counts)
# 8546  211

#### End of optional: Dropout removal




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
table (seurat.ss2@meta.data$label)   # run_de will compare G2 vs G4
# G2  G4 
#109 102 

mymean <- data.frame (mean= apply (counts, 1, mean))


res <- run_de(seurat.ss2, de_method = 'wilcox', de_family= "singlecell")
res <- data.frame (res)
row.names (res) <- res$gene
res <- merge (res, mymean, by="row.names")
res <- res[order (res$p_val_adj), ]
row.names (res) <- res$gene
colnames (res)[1] <- "gene_name"
res$avg_logFC <- -1 * res$avg_logFC


## Add the cell proportion
res <- merge (res,prop, by="row.names")
res <- res[ ,-1]
res <- res[ ,-dim (res)[2]]
row.names (res) <- res$gene_name 


## Annotation

library('org.Mm.eg.db')

#columns(org.Mm.eg.db)
symbols <- row.names (res)
res1a <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (row.names (res), names (res1a))
res$Description <- as.vector (res1a) [idx]
res <- res[order (res$p_val_adj), ]
res <- res[ ,-which (colnames (res) == "de_family")]

write.xlsx (res, "table 3. hippocampus_G4vsG2_selected_cells_normalization_wilcoxon_analysis.xlsx", rowNames=F)











