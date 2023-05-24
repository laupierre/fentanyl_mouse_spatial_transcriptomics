## hippocampus of G1 (Sham+Veh) vs G2 (SNI+Veh)  

# devtools::install_github("neurorestore/Libra")
library (Libra)
library (Seurat)
library (openxlsx)

options(Seurat.object.assay.version = "v5")


## Wilcoxon test
## Here we are normalizing on the selected cells only after extracting the raw counts from these cells

# str (brain@assays)
# raw counts
# counts <- data.matrix (brain[["Spatial"]]$counts)
# max (as.data.frame (counts [ ,1]))

# normalized counts (log2)
# counts <- as.matrix (brain@assays$SCT@data) 
# max (as.data.frame (counts [ ,1]))

# scaled data
# counts <- as.matrix (brain@assays$SCT@scale.data) 
# max (as.data.frame (counts [ ,1]))


brain <- readRDS ("brain_G2G1_groups.rds")

# raw counts
counts <- as.matrix (brain[["SCT"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])
dim (counts)
# 18827   229

boxplot (apply (counts, 1, mean))
counts <- counts[apply (counts, 1, mean) > 0.05, ]
dim (counts)
# 12381   229



#### Optional: Dropout analysis

meta <- brain@meta.data
meta.g1 <- meta[grep ("G1", meta$group), ]
meta.g1 <- meta.g1 [row.names (meta.g1) %in% colnames (counts), ]
table (meta.g1$group)

meta.g2 <- meta[grep ("G2", meta$group), ]
meta.g2 <- meta.g2 [row.names (meta.g2) %in% colnames (counts), ]
table (meta.g2$group)


counts.g1 <- counts[ ,colnames (counts) %in% row.names (meta.g1)]
dim (counts.g1)
counts.g2 <- counts[ ,colnames (counts) %in% row.names (meta.g2)]
dim (counts.g2)


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
  ggtitle ("Percentage of zeros in the G1 group")
  
p2 <- ggplot(data = prop, mapping = aes(x=prop2)) + 
  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7, bins=20)  +
  geom_density() + xlab ("Percentage of zeros") + ylab ("Density of expressed genes") + geom_vline(xintercept=75, linetype="dashed", color = "red") +
  ggtitle ("Percentage of zeros in the G3 group") 
  
p3 <- ggarrange (p1, p2, nrow=1)
p3


# Proportion of cells expressing a gene
prop <- 100- prop

# select genes with at least 25% of cells expressing the gene in any group
prop$drop <- apply (prop, 1, function (x) {any (x > 25)})
prop <- prop[prop$drop == TRUE, ]

counts <- counts[row.names (counts) %in% row.names (prop), ]
dim (counts)
# 9199  229

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
table (seurat.ss2@meta.data$label)   # It will compare G1 vs G2
# G1  G2 
#120 109 

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

write.xlsx (res, "table 1. hippocampus_G2vsG1_selected_cells_normalization_wilcoxon_analysis.xlsx", rowNames=F)

boxplot (res$avg_logFC)
abline (h=0)



#########
## Verification of the log fold change orientation using voom-limma

library (edgeR)
library (limma)

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
  
condition <- meta.ss2$label
mm <- model.matrix(~0 + condition)
y <- voom(d0, mm, plot = F)
fit <- lmFit(y, mm)

# G2 (SNI+Veh) vs G1 (Sham+Veh)
contr <- makeContrasts(conditionG2 - conditionG1, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res.lim <- topTable(tmp, sort.by = "p", n = Inf) 
res.lim <- res.lim[res.lim$adj.P.Val <= 0.05, ]

resa <- merge (res, res.lim, by="row.names")


## With this G2 vs G1 comparison, we have the same trend
par (mfrow=c(2,1))
plot (resa$avg_logFC, resa$logFC, xlab="Wilcoxon log fold change", ylab="limma log fold change (G2 vs G1)")
abline (h=0)
abline (v=0)

plot (resa$avg_logFC, -resa$logFC, xlab="Wilcoxon log fold change", ylab="limma log fold change (G1 vs G2)")
abline (h=0)
abline (v=0)
























