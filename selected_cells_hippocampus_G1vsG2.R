## hippocampus of G1 (Sham+Veh) vs G2 (SNI+Veh)  

library (Seurat)

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

brain <- merge(brain1, y = c(brain2, brain3, brain4), add.cell.ids = c("2C", "2A", "1C", "1A"))
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
saveRDS (brain, "brain_G2G1_groups.rds")


#counts <- as.matrix (brain1@assays$SCT@data [ ,WhichCells(brain1, expression = location == "Hippocampus")])
#dim (counts)














