## hippocampus of G1 (Sham+Veh) vs G2 (SNI+Veh)  

library (Seurat)
library (openxlsx)

options(Seurat.object.assay.version = "v5")

# generic function
preprocess <- function (data.dir, meta) {
sample.name <- gsub (".*/", "", data.dir)
brain <- Load10X_Spatial (data.dir, filename= "filtered_feature_bc_matrix.h5")
brain[["Spatial"]] <- as(brain[["Spatial"]], Class = "Assay5")

# Remove cells with mitochondrial contamination
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "percent_mito")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")
brain <- brain[ ,brain$nFeature_Spatial > 500 & brain$percent_mito < 20 & brain$percent_hb < 20]

# remove mitochondrial genes
brain <- brain[!grepl("^mt-", rownames(brain)), ]
# remove Hemoglobin genes (if that is a problem on your data)
brain <- brain[!grepl("^Hb.*-", rownames(brain)), ]

# normalization
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

colnames (meta)[2] <- "TOTAL"
meta <- meta[meta$TOTAL == "rHEMI", ]
meta1 <- meta[meta$Hip != "", ]
print (dim (meta1))

cells.use <- meta1$Barcode
cells.use <- cells.use[cells.use %in% row.names (brain@meta.data)]
idx <- match (cells.use, row.names (brain@meta.data))
brain@meta.data$location <- "unknown"
brain@meta.data$location [idx] <- "Hippocampus"
brain@meta.data$group <- sample.name
brain@meta.data$cell <- paste (row.names(brain@meta.data), sample.name, sep="-")

brain <- UpdateSeuratObject(brain)
return (brain)
}


data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G2-2C"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G2_2C_sniv02.csv", sep=",")
brain1 <- preprocess (data.dir, meta)

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G2-2A"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G2_2A_sniv01.csv", sep=",")
brain2 <- preprocess (data.dir, meta)

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1C"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1C_shmv02.csv", sep=",")
brain3 <- preprocess (data.dir, meta)

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01.csv", sep=",")
brain4 <- preprocess (data.dir, meta)


# By default, merge() will combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices. 
# If you want to merge the normalized data matrices as well as the raw count matrices, simply pass merge.data = TRUE. 

brain <- merge(brain1, y = c(brain2, brain3, brain4), add.cell.ids = c("2C", "2A", "1C", "1A"))

DefaultAssay(brain) <- "Spatial"
brain <- JoinLayers (brain)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
saveRDS (brain, "brain_G2G1_groups.rds")



################
## Wilcoxon test
## Here we are normalizing on the selected cells only after extracting the raw counts from these cells

# devtools::install_github("neurorestore/Libra")
library (Libra)

# str (brain@assays)
# raw counts
# counts <- data.matrix (brain[["Spatial"]]$counts)
# max (as.data.frame (counts [ ,1]))

# normalized counts
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


## Annotation

library('org.Mm.eg.db')

#columns(org.Mm.eg.db)
symbols <- row.names (res)
res1a <- mapIds(org.Mm.eg.db, symbols, 'GENENAME', 'SYMBOL')

idx <- match (row.names (res), names (res1a))
res$Description <- as.vector (res1a) [idx]
res <- res[order (res$p_val_adj), ]
res <- res[ ,-which (colnames (res) == "de_family")]

write.xlsx (res, "hippocampus_G2vsG1_selected_cells_wilcoxon_analysis.xlsx", rowNames=F)

boxplot (res$avg_logFC)
abline (h=0)


#########

brain.n <- read.xlsx ("hippocampus_G2vsG1_selected_cells_brain_normalization_wilcoxon_analysis.xlsx")

comp2 <- merge (res, brain.n, by="gene_name") 

pdf ("Comparison hippocampus and normalization methods G2 vs G1.pdf")
plot (comp2$log.fold.change, comp2$avg_logFC, xlab="log fold changes wilcoxon (brain norm)", ylab="log fold changes wilcoxon (hippocampus norm)", main="Comparison hippocampus and normalization methods",
      xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj < 0.05, "blue","black"))
abline (0,1, col="red")
abline (h=0)
abline (v=0)
dev.off ()



#########
## Sanity check (with bulk RNA-Seq)

rnaseq <- read.xlsx ("sham_vs_sni_Differential_Expression.xlsx")
rnaseq <- rnaseq[ ,c("Gene.Symbol", "Log2.Fold.Change", "FDR.Adj.p.Value", "Mean")]

res2 <- merge (res, rnaseq, by.x="gene_name", by.y="Gene.Symbol")
head (res2)

pdf ("Comparison spatial and RNA-Seq methods G2 vs G1.pdf")
plot (res2$avg_logFC, res2$Log2.Fold.Change, xlab="Spatial log fold changes Wilcoxon", ylab="Chicago log fold changes RNA-Seq", main="Comparison spatial and RNA-Seq methods", 
     xlim=c(-2,2), ylim=c(-2,2), col=ifelse (comp2$padj < 0.05, "blue","black"))
abline (h=0)
abline (v=0)
dev.off ()



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
























