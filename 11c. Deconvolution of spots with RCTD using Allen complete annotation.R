library (openxlsx)
library (spacexr)
library (Seurat)
library (dplyr)

# load the visium data
dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A"
visium <- read.VisiumSpatialRNA (dir)

# Create a list of barcodes from the column names of the count matrix
barcodes <- colnames(visium@counts)

# Plot number of UMIs per barcode (spot). nUMI is the nCount_Spatial (initially present in the Seurat metadata)
## The number of UMIs (RNA species) is proportional to the number of cells 

plot_puck_continuous(puck=visium, barcodes=barcodes, plot_val=visium@nUMI, 
                     size=1, ylimit=c(0,round(quantile(visium@nUMI,0.9))), 
                     title='plot of nUMI') 
       


## see allen_mouse_adult_brain_single_cell repository (Github) for Allen 1 million cells reference

allen_reference <- readRDS("seurat_allen_small_mouse_brain_with_UMAP.RDS")

Idents (allen_reference) <- "subclass_label"
DimPlot(allen_reference, label = TRUE)





