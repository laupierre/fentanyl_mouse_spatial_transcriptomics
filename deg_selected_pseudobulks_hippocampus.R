# /Volumes/king/haram2023_spatial

library (Libra)
library (Seurat)
library (edgeR)
library (limma)
library (openxlsx)

brain <- readRDS ("brain_G2G1_groups.rds")

# raw counts
counts <- as.matrix (brain[["Spatial"]]$counts [ ,WhichCells(brain, expression = location == "Hippocampus")])
dim (counts)
# 32264   229

boxplot (apply (counts, 1, mean))
counts <- counts[apply (counts, 1, mean) > 0.5, ]
dim (counts)
# 8061  229


meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == "Hippocampus"), ]
# we need cell_type (hippocampus), replicate (mouse origin), and label (treatment, i.e G2 or G1 groups)
meta.ss2$replicate <- gsub (".*-", "", meta.ss2$group)
meta.ss2$label <- gsub ("-.*", "", meta.ss2$group)
meta.ss2$cell_type <- "Hippocampus"
meta <- meta.ss2

idx <- match (colnames (counts), row.names (meta))
meta <- meta[idx, ]
stopifnot (colnames (counts) == row.names (meta))

mymean <- data.frame (mean= apply (counts, 1, mean))



##  Pseudobulk construction

meta$mouse <- paste (meta$replicate, meta$label, sep=":")
mouse <- unique (meta$mouse)
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



## Limma chunk (normal voom)

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
# 5
res.voom <- res

res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- cbind (data.frame (gene_name= row.names (res), res))

write.xlsx (res, "pseudobulk_hippocampus_selected_cells.xlsx", rowNames=F)



## No main differences when compared to the alternative methods 
## voomQW with sample variability
y <- voomWithQualityWeights(d0, design = mm, plot = TRUE)
fit <- lmFit(y, mm)
contr <- makeContrasts(conditionG2 - conditionG1 , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
res <- topTable(tmp, sort.by = "p", n = Inf) 
res <- res[res$adj.P.Val < 0.05, ]

table (res$adj.P.Val < 0.05)
# 5
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
# 5
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
# 4
res.voomG <- res






