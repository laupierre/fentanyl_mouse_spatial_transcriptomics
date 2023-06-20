# /Volumes/king/haram2023_spatial

library (Libra)
library (Seurat)
library (edgeR)
library (limma)
library (openxlsx)
library (sctransform)

brain <- readRDS ("brain_G2G1_groups.rds")


## Here we are normalizing the hippocampus cells only with lognorm (and not the entire brain)

#counts <- as.matrix (brain[["SCT"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])

myarea <- "DG"
counts <- as.matrix (brain[["Spatial"]]$counts [ ,WhichCells(brain, expression = location == myarea)])
dim (counts)
# 18827   229

# take the corrected counts (v2)
# brain <- SCTransform(brain, vst.flavor = "v2") 
# counts <- as.matrix (brain[["SCT"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])
# dim (counts)
# 18827   229

boxplot (apply (counts, 1, mean))
counts <- counts[apply (counts, 1, mean) > 0.05, ]
dim (counts)
# 12381   229


#meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == "Hippocampus"), ]
meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == myarea), ]
# we need cell_type (hippocampus), replicate (mouse origin), and label (treatment, i.e G2 or G1 groups)
meta.ss2$replicate <- gsub (".*-", "", meta.ss2$group)
meta.ss2$label <- gsub ("-.*", "", meta.ss2$group)
#meta.ss2$cell_type <- "Hippocampus"
meta.ss2$cell_type <- myarea
meta <- meta.ss2

idx <- match (colnames (counts), row.names (meta))
meta <- meta[idx, ]
stopifnot (colnames (counts) == row.names (meta))

mymean <- data.frame (mean= apply (counts, 1, mean))



##  Pseudobulk construction of 4 pseudo-mice

meta$mouse <- paste (meta$replicate, meta$label, sep=":")
mouse <- unique (meta$mouse)
mouse
# "2C:G2" "2A:G2" "1C:G1" "1A:G1"


mat <- list ()
for (i in (1:length (mouse))) {
mycells <- counts[ ,colnames (counts) %in% row.names (meta)[meta$mouse == mouse[i] ]]
a <- data.frame (rowSums (mycells))
row.names (a) <- row.names (mycells)
colnames (a)[1] <- mouse[i]
mat[[i]] <- a
}

pseudo.counts <- do.call ("cbind", mat)
head (pseudo.counts)



## Limma chunk (use the normal voom to normalize the hippocampal pseudobulks) 

d0 <- DGEList(pseudo.counts)
d0 <- calcNormFactors(d0)

condition <- gsub (".*:", "", colnames (pseudo.counts))
mm <- model.matrix(~0 + condition)

y <- voom(d0, mm, plot = T)

fit <- lmFit(y, mm)
contr <- makeContrasts(conditionG2 - conditionG1 , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- res[res$adj.P.Val <= 0.05, ]

table (res$adj.P.Val < 0.05)
# 2 vs 44
res.voom <- res

res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- cbind (data.frame (gene_name= row.names (res), res))

write.xlsx (res, paste (myarea, "_area_pseudobulk_hippocampus_lognorm_selected_cells.xlsx", sep=""), rowNames=F)


## Sanity check

sel <- read.xlsx (paste (paste ("table 1.", myarea), "_area_G2vsG1_selected_cells_normalization_wilcoxon_analysis_with_percentage_cells.xlsx", sep=""))

sel <- merge (sel, res, by.x="gene", by.y="gene_name")
head (sel)

plot (sel$avg_logFC, sel$logFC, xlab= "normalization selected cells (Libra)", ylab="Pseudobulks")
abline (0,1,col="red")
abline (h=0)
abline (v=0)




############
##### No main differences when compared to the alternative methods 

## voomQW with sample variability
y <- voomWithQualityWeights(d0, design = mm, plot = TRUE)
fit <- lmFit(y, mm)
contr <- makeContrasts(conditionG2 - conditionG1 , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- res[res$adj.P.Val < 0.05, ]

table (res$adj.P.Val < 0.05)
# 2 vs 43
res.voomQW <- res



## voomQW with block variability
y <- voomWithQualityWeights(d0, design = mm, var.group=condition , plot = TRUE)
fit <- lmFit(y, mm)
contr <- makeContrasts(conditionG2 - conditionG1 , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- res[res$adj.P.Val < 0.05, ]

table (res$adj.P.Val < 0.05)
# 2 vs 43
res.voomQWB <- res



## voomQW with group variability
## see the group-specific mean-variance plots to assess the usage of this method
# system ("git clone https://github.com/YOU-k/voomByGroup.git")
source ("~/voomByGroup/voomByGroup.R")

y <- voomByGroup(d0, design = mm, group=condition , plot = "all")
fit <- lmFit(y, mm)
contr <- makeContrasts(conditionG2 - conditionG1 , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- res[res$adj.P.Val < 0.05, ]

table (res$adj.P.Val < 0.05)
# 2 vs 32
res.voomG <- res






