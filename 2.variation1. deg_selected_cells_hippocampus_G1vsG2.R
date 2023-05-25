## hippocampus of G1 (Sham+Veh) vs G2 (SNI+Veh)  

# devtools::install_github("neurorestore/Libra")
library (Libra)
library (Seurat)
library (openxlsx)

options(Seurat.object.assay.version = "v5")


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

## Step1 (numa)

meta <- brain@meta.data
meta.g1 <- meta[grep ("G1-1A", meta$group), ]
meta.g1 <- meta.g1 [row.names (meta.g1) %in% colnames (counts), ]
table (meta.g1$group)

meta.g2 <- meta[grep ("G2-2A", meta$group), ]
meta.g2 <- meta.g2 [row.names (meta.g2) %in% colnames (counts), ]
table (meta.g2$group)

meta.g3 <- meta[grep ("G1-1C", meta$group), ]
meta.g3 <- meta.g3 [row.names (meta.g3) %in% colnames (counts), ]
table (meta.g3$group)

meta.g4 <- meta[grep ("G2-2C", meta$group), ]
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

# Percentage of cells expressing a gene
prop1 <- data.frame (G1_1A= apply (counts.g1, 1, function (x) {sum (x == 0)}))
prop2 <- data.frame (G2_2A= apply (counts.g2, 1, function (x) {sum (x == 0)}))
prop3 <- data.frame (G1_1C= apply (counts.g3, 1, function (x) {sum (x == 0)}))
prop4 <- data.frame (G2_2C= apply (counts.g4, 1, function (x) {sum (x == 0)}))
prop1  <- table (meta.g1$group)[[1]] - prop1
prop2  <- table (meta.g2$group)[[1]] - prop2
prop3  <- table (meta.g3$group)[[1]] - prop3
prop4  <- table (meta.g4$group)[[1]] - prop4
numa <- cbind (prop1, prop2, prop3, prop4)
numa <- numa[ ,order (colnames (numa))]

prop1a  <- (prop1 / table (meta.g1$group)[[1]] ) *100
prop2a  <- (prop2 / table (meta.g2$group)[[1]] ) *100
prop3a  <- (prop3 / table (meta.g3$group)[[1]] ) *100
prop4a  <- (prop4 / table (meta.g4$group)[[1]] ) *100
numa_a <- round (cbind (prop1a, prop2a, prop3a, prop4a))
colnames (numa_a) <- paste (colnames (numa_a), ".pct", sep="")
numa_a <- numa_a[ ,order (colnames (numa_a))]

numa <- cbind (numa, numa_a)
head (numa)


# Compute two-proportions z-test
# Fisher Exact probability test is an excellent non-parametric technique for comparing proportions, when the two independent samples are small in size.
# We want to know, whether the proportions of positive spots are the same in the two groups (G1 and G2)
# See http://www.sthda.com/english/wiki/two-proportions-z-test-in-r

resl <- list ()

for (i in (1:dim (numa)[1])) {
x <- c (numa[i, 1] + numa[i ,2], numa[i, 3] + numa[i, 4])
n <- c (table (meta.g1$group)[[1]] + table (meta.g3$group)[[1]] , table (meta.g2$group)[[1]] + table (meta.g4$group)[[1]] )

res <- prop.test(x = x, n = n)
res <- cbind (t (data.frame (res$estimate)) *100, data.frame (pvalz= res$p.value))
row.names (res) <- row.names (numa)[i]
resl[[i]] <- res
}

resl <- do.call ("rbind", resl)
resl$pvalz_adj <- p.adjust (resl$pvalz, method= "BH")
head (resl)












