library (Seurat)

options(Seurat.object.assay.version = "v5")

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A"
#meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G2_2C_sniv02.csv", sep=",")
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A/cell_proportions_allen.txt")


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
#brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

colnames (meta)[1] <- "Barcode"
meta1 <- meta
print (dim (meta1))

tokeep <- row.names (brain@meta.data) [row.names (brain@meta.data) %in% meta1$Barcode]
brain <- subset (brain, cells = tokeep)
idx <- match (tokeep, row.names (brain@meta.data))
brain@meta.data$allen <- meta1$cell_type1[idx]
brain@meta.data$allen <- gsub ("=.*","", brain@meta.data$allen)
brain@meta.data$group <- sample.name
brain@meta.data$cell <- paste (row.names(brain@meta.data), sample.name, sep="-")

brain <- UpdateSeuratObject(brain)
return (brain)

Idents (brain) <- "allen"
SpatialDimPlot(brain)
#SpatialDimPlot(brain, group.by = c("allen"))

dg <- subset(brain, idents = c("DG"))
SpatialDimPlot(dg)









