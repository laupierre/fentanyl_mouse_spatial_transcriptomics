## Slide 1 containing the G1 and G3 samples

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
#brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

colnames (meta)[2] <- "TOTAL"
meta <- meta[meta$TOTAL == "rHEMI", ]
meta1 <- meta[meta$Hip != "", ]
print (dim (meta1))

cells.use <- meta1$Barcode
cells.use <- cells.use[cells.use %in% row.names (brain@meta.data)]
#idx <- match (cells.use, row.names (brain@meta.data))
# new
idx <- match (meta1$Barcode, row.names (brain@meta.data))
brain@meta.data$location <- "unknown"
#brain@meta.data$location [idx] <- "Hippocampus"
# new
brain@meta.data$location[idx] <- meta1$Hip   
brain@meta.data$group <- sample.name
brain@meta.data$cell <- paste (row.names(brain@meta.data), sample.name, sep="-")
# sanity check
san <- merge (brain@meta.data, meta1, by.x="row.names", by.y="Barcode")
stopifnot (san$location == san$Hip)
brain <- UpdateSeuratObject(brain)
  
return (brain)
}


data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01.csv", sep=",")
brain1 <- preprocess (data.dir, meta)

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G3-1B"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G3_1B_shmf01.csv", sep=",")
brain2 <- preprocess (data.dir, meta)

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1C"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1C_shmv02.csv", sep=",")
brain3 <- preprocess (data.dir, meta)

data.dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G3-1D"
meta <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G3_1D_shmf02.csv", sep=",")
brain4 <- preprocess (data.dir, meta)


# By default, merge() will combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices. 
# If you want to merge the normalized data matrices as well as the raw count matrices, simply pass merge.data = TRUE. 

brain <- merge(brain1, y = c(brain2, brain3, brain4), add.cell.ids = c("1A", "1B", "1C", "1D"))

DefaultAssay(brain) <- "Spatial"
brain <- JoinLayers (brain)
## SCT normalization of all slides (deprecated method)
#brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
saveRDS (brain, "brain_slide1_G1G3_groups.rds")

