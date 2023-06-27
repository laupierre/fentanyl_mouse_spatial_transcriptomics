## Differential gene expression with option 2. deg_selected_cells_hippocampus_G1vsG2.R 


## hippocampus of G1 (Sham+Veh) vs G2 (SNI+Veh)  

# devtools::install_github("neurorestore/Libra")
library (Libra)
library (Seurat)
library (openxlsx)

options(Seurat.object.assay.version = "v5")

brain <- readRDS ("brain_G2G1_groups.rds")

myarea <- "DG"
counts <- as.matrix (brain[["Spatial"]]$counts [ ,WhichCells(brain, expression = location == myarea)])
dim (counts)
# 18827   229

boxplot (apply (counts, 1, mean))
counts <- counts[apply (counts, 1, mean) > 0.05, ]
dim (counts)


## Step 2: Gene dropout removal 
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
  ggtitle ("Percentage of zeros in the G2 group") 
  
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

#meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == "Hippocampus"), ]
meta.ss2 <- brain@meta.data[WhichCells(brain, expression = location == myarea), ]
meta.ss2$group.use <- gsub ("-.*", "", meta.ss2$group)

seurat.ss2@meta.data <- meta.ss2
head (seurat.ss2@meta.data)


# Here we are normalizing to populate the data slot !!!
seurat.ss2 <- NormalizeData (seurat.ss2)


library(presto)

# here, we a comparing G2 and G1 groups, however, both directions are reported !
res_AB <- wilcoxauc (seurat.ss2, "group.use", assay = 'data', groups_use = c('G2', 'G1'))
res_AB <- res_AB[order (res_AB$padj), ]
res_AB <- res_AB[res_AB$group == "G2", ]
head (res_AB)



#### fgsea

library(msigdbr)
library(fgsea)
library(dplyr)
library (tidyverse)
library(ggplot2)
library (openxlsx)

msigdbr_species()

m_df <- msigdbr(species = "Mus musculus", category = "C5")
head(m_df)

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head (fgsea_sets)


## select only the feature and auc columns for fgsea, which statistics to use is an open question !

genes <- res_AB %>% dplyr::filter(group == "G2") %>%
                    arrange(desc(auc)) %>% 
                    dplyr::select(feature, auc)

head (genes)

ranks <- deframe(genes)
head(ranks)

fgseaRes <- fgseaMultilevel (fgsea_sets, stats = ranks, minSize=15, maxSize=500)

fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
head (data.frame (fgseaResTidy))

write.xlsx (data.frame (fgseaResTidy), paste (myarea, "_G2vsG1_fgsea_results.xlsx"), rowNames=F)


p <- plotEnrichment(fgsea_sets[["GOCC_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX"]], ranks) + ggtitle ("GOCC_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX")
ggsave ("fgsea_plot.pdf", p)










