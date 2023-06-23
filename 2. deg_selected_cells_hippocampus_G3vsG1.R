library (Libra)
library (Seurat)
library (openxlsx)

options(Seurat.object.assay.version = "v5")

brain <- readRDS ("brain_slide1_G1G3_groups.rds")

# raw counts
#counts <- as.matrix (brain[["SCT"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])
#dim (counts)
# 18966   252

myarea <- "DG"
counts <- as.matrix (brain[["Spatial"]]$counts [ ,WhichCells(brain, expression = location == myarea)])
dim (counts)

counts <- counts[apply (counts, 1, mean) > 0.05, ]
dim (counts)
# 12691   252



#### Optional: Dropout removal

## Step1

meta <- brain@meta.data

meta.g1 <- meta[grep ("G1-1A", meta$group), ]
meta.g1 <- meta.g1 [row.names (meta.g1) %in% colnames (counts), ]
table (meta.g1$group)

meta.g2 <- meta[grep ("G3-1B", meta$group), ]
meta.g2 <- meta.g2 [row.names (meta.g2) %in% colnames (counts), ]
table (meta.g2$group)

meta.g3 <- meta[grep ("G1-1C", meta$group), ]
meta.g3 <- meta.g3 [row.names (meta.g3) %in% colnames (counts), ]
table (meta.g3$group)

meta.g4 <- meta[grep ("G3-1D", meta$group), ]
meta.g4 <- meta.g4 [row.names (meta.g4) %in% colnames (counts), ]
table (meta.g4$group)


counts.g1 <- counts[ ,colnames (counts) %in% row.names (meta.g1)]
dim (counts.g1)
counts.g2 <- counts[ ,colnames (counts) %in% row.names (meta.g2)]
dim (counts.g2)
counts.g3 <- counts[ ,colnames (counts) %in% row.names (meta.g3)]
dim (counts.g3)
counts.g4 <- counts[ ,colnames (counts) %in% row.names (meta.g4)]
dim (counts.g4)

# Number of cells that express a gene (for all genes)
prop1 <- data.frame (G1_1A= apply (counts.g1, 1, function (x) {sum (x != 0)}))
prop2 <- data.frame (G3_1B= apply (counts.g2, 1, function (x) {sum (x != 0)}))
prop3 <- data.frame (G1_1C= apply (counts.g3, 1, function (x) {sum (x != 0)}))
prop4 <- data.frame (G3_1D= apply (counts.g4, 1, function (x) {sum (x != 0)}))
numa <- cbind (prop1, prop2, prop3, prop4)

# Add the two numbers of cells
prop <- cbind ( data.frame (G1= prop1+prop3), data.frame (G3= prop2+prop4))
colnames (prop) <- c("G1","G3")

numa <- cbind (numa, prop)

# Percentage of expressing cells in a group
prop2 <- data.frame (G1_prop= prop$G1 / (dim (counts.g1) [2] + dim (counts.g3) [2]))
prop3 <- data.frame (G3_prop= prop$G3 / (dim (counts.g2) [2] + dim (counts.g4) [2]))
prop <- round (cbind (prop2 *100, prop3*100), digits=1)

numa <- cbind (numa, prop)
head (numa)



## Step2 
meta <- brain@meta.data
meta.g1 <- meta[grep ("G1", meta$group), ]
meta.g1 <- meta.g1 [row.names (meta.g1) %in% colnames (counts), ]
table (meta.g1$group)

meta.g2 <- meta[grep ("G3", meta$group), ]
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

#### End of optional step: Gene dropout removal




# for backward compatibility
counts <- as (counts, "sparseMatrix")
seurat.ss2 <- CreateSeuratObject(counts)

meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == myarea), ]
# we need cell_type (eg hippocampus), replicate (eg mouse), and label (eg treatment)
meta.ss2$replicate <- gsub (".*-", "", meta.ss2$group)
meta.ss2$label <- gsub ("-.*", "", meta.ss2$group)
meta.ss2$cell_type <- myarea

meta.ss2 <- meta.ss2[row.names (meta.ss2) %in% colnames (counts), ]
idx <- match (row.names (meta.ss2), colnames (counts))
meta.ss2 <- meta.ss2[idx, ]
stopifnot (row.names (meta.ss2) == colnames (counts))

seurat.ss2@meta.data <- meta.ss2
table (seurat.ss2@meta.data$label)   # run_de will compare G1 vs G3
# G1  G3 
#120 132 


res <- run_de(seurat.ss2, de_method = 'wilcox', de_family= "singlecell")
res <- data.frame (res)
row.names (res) <- res$gene
# we want G3 vs G1
res$avg_logFC <- -1 * res$avg_logFC

## Add numa info
res <- merge (res, numa, by.x="gene", by.y="row.names")
row.names (res) <- res$gene 
res <- res[order (res$p_val_adj), ]

## Annotation
library('org.Mm.eg.db')

symbols <- row.names (res)
res1a <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (row.names (res), names (res1a))
res$Description <- as.vector (res1a) [idx]
res <- res[order (res$p_val_adj), ]
res <- res[ ,-which (colnames (res) == "de_family")]

table (res$p_val_adj < 0.05)
#FALSE  TRUE 
# 9963     3 

boxplot (res$avg_logFC)
abline (h=0)

write.xlsx (res, paste (paste ("table 2.", myarea), "_area_G3vsG1_selected_cells_normalization_wilcoxon_analysis_with_percentage_cells.xlsx", sep=""), rowNames=F)






