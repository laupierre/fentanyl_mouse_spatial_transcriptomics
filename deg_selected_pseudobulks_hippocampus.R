# /Volumes/king/haram2023_spatial

library (Libra)
library (Seurat)
library (edgeR)
library (limma)


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
# we need cell_type, replicate, and label
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







library (edgeR)
library (limma)

d0 <- DGEList(pseudo.counts)
d0 <- calcNormFactors(d0)
  
condition <- c (rep ("A", 3), rep ("B", 3))

mm <- model.matrix(~0 + condition)

y <- voom(d0, mm, plot = T)


