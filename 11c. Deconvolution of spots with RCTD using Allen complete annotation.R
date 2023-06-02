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


allen_reference@active.assay = "RNA"

# wilcoxon test between annotated clusters
markers_sc <- FindAllMarkers(allen_reference, only.pos = TRUE, logfc.threshold = 0.1,
              test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
              return.thresh = 0.05, assay = "RNA")

# Select genes that are also present in the ST data (in addition to the Allen SC data)
markers_sc <- markers_sc[markers_sc$gene %in% rownames(visium@counts), ]

# Select top 20 genes per cluster, select top by first p-value, then absolute diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 + 1))
markers_sc %>%
    group_by(cluster) %>%
    top_n(-100, p_val) %>%
    top_n(50, pct.diff) %>%
    top_n(20, log.pct.diff) -> top20

m_feats <- unique(as.character(top20$gene))


# Create raw counts
allen_reference <- allen_reference[row.names (allen_reference) %in% m_feats, ]
counts <- allen_reference[["RNA"]]@data
dim (counts)
# 579 8152

# See https://www.10xgenomics.com/resources/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
# Create factors of cell types
cell_types_idx <- match (colnames (counts), names (allen_reference$subclass_label))
cell_types <- allen_reference$subclass_label [cell_types_idx]

# Remove / in names
cell_types[cell_types == "L2/3 IT CTX"] <- "L2_3 IT"
cell_types[cell_types == "L2/3 IT ENTl"] <- "L2_3 IT ENTl"
cell_types[cell_types == "L2/3 IT PPP"] <- "L2_3 IT PPP"
cell_types[cell_types == "L2/3 IT RHP"] <- "L2_3 IT RHP"
cell_types[cell_types == "L4/5 IT CTX"] <- "L4_5 IT CTX"
cell_types[cell_types == "L5/6 IT TPE-ENT"] <- "L5_6 IT TPE-ENT"
cell_types[cell_types == "L5/6 NP CTX"] <- "L5_6 NP CTX"
cell_types[cell_types == "L6b/CT ENT"] <- "L6b_CT ENT"

cell_types <- factor (cell_types)
table (cell_types)


# Create nUMI
nUMI <- allen_reference@meta.data$nCount_RNA
nUMI <- nUMI[cell_types_idx]
names (nUMI) <- colnames (counts)

reference <- Reference(counts, cell_types, nUMI)
# Reference: some nUMI values are less than min_UMI = 100, and these cells will be removed


## RCTD in multimode

Sys.setenv("OPENBLAS_NUM_THREADS"=2)
myRCTD <- create.RCTD(visium, reference, max_cores = 2)

# full mode, which assigns any number of cell types per spot and is recommended for technologies with poor spatial resolution such as 100-micron resolution Visium; 
# multi mode, an extension of doublet mode that can discover more than two cell types per spot as an alternative option to full mode

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
saveRDS (myRCTD, "myRCTD_visium_allen_full_annot_multi_mode.rds")








