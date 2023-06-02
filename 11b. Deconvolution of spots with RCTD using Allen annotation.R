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
           
           

## get the Allen brain mouse cortex annotation
## see https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Select_genes_for_deconvolution

## Keep in mind that it is important to have a reference that contains all the celltypes you expect to find in your spots. 
## Ideally it should be a scRNAseq reference from the exact same tissue.
## We will use a reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol.


# system ("wget -O allen_cortex.rds https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1")

allen_reference <- readRDS("allen_cortex.rds")

# check number of cells per subclass (there is no hippocampus there !!!)
table(allen_reference$subclass)

# select 200 cells max per subclass, first set subclass as active.ident
Idents(allen_reference) <- allen_reference$subclass
allen_reference <- subset(allen_reference, cells = WhichCells(allen_reference, downsample = 200))

# check again number of cells per subclass
table(allen_reference$subclass)

# run SCT normalization and dimensionality reduction
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE, method = "poisson") %>%
                   RunPCA(verbose = FALSE) %>%
                   RunUMAP(dims = 1:30)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, label = TRUE)



