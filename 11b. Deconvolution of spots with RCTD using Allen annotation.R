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
           
           

## get the standard Allen brain mouse cortex annotation
## see https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Select_genes_for_deconvolution

## Keep in mind that it is important to have a reference that contains all the celltypes you expect to find in your spots. 
## Ideally it should be a scRNAseq reference from the exact same tissue.
## We will use a reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol.


# system ("wget -O allen_cortex.rds https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1")

allen_reference <- readRDS("allen_cortex.rds")

# check number of cells per subclass (there is no hippocampus there !!!)
table(allen_reference$subclass)

# select 100 cells max per subclass, first set subclass as active.ident
Idents(allen_reference) <- allen_reference$subclass
allen_reference <- subset(allen_reference, cells = WhichCells(allen_reference, downsample = 100))

# check again number of cells per subclass
table(allen_reference$subclass)

# remove cell types with count lower than 25 (CR neurons)
torm <- WhichCells(allen_reference, ident = "CR")
allen_reference <- allen_reference[ ,!colnames (allen_reference) %in% torm]
table(allen_reference$subclass)

# run SCT normalization and dimensionality reduction
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE, method = "poisson") %>%
                   RunPCA(verbose = FALSE) %>%
                   RunUMAP(dims = 1:30)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, label = TRUE)


# Deconvolution is a method to estimate the abundance (or proportion) of different celltypes in a bulkRNAseq dataset using a single cell reference.
# Most deconvolution methods does a prior gene selection and there are different options that are used:
# Use variable genes in the SC data.
# Use variable genes in both SC and ST data
# DE genes between clusters in the SC data. (what we do below)

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
# 384 1862

# See https://www.10xgenomics.com/resources/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
# Create factors of cell types

cell_types_idx <- match (colnames (counts), names (allen_reference$subclass))
cell_types <- allen_reference$subclass [cell_types_idx]

# Remove / in names
cell_types[cell_types == "L2/3 IT"] <- "L2_3 IT"

cell_types <- factor (cell_types)
table (cell_types)


# Create nUMI
nUMI <- allen_reference@meta.data$nCount_RNA
nUMI <- nUMI[cell_types_idx]
names (nUMI) <- colnames (counts)

reference <- Reference(counts, cell_types, nUMI)


## RCTD in multimode

Sys.setenv("OPENBLAS_NUM_THREADS"=2)
myRCTD <- create.RCTD(visium, reference, max_cores = 2)

# full mode, which assigns any number of cell types per spot and is recommended for technologies with poor spatial resolution such as 100-micron resolution Visium; 
# multi mode, an extension of doublet mode that can discover more than two cell types per spot as an alternative option to full mode

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
saveRDS (myRCTD, "myRCTD_visium_allen_multi_mode.rds")


# The results of RCTD multimode are stored in myRCTD@results
# Each entry represents the estimated proportion of each cell type on each pixel.

propl <- list ()

for (i in (1:length (myRCTD@results))) {
pixel <- myRCTD@results[[i]]
prop <- t (data.frame (pixel$all_weights))
ratio <- 1/sum (prop) 
prop <- round (prop * ratio * 100)
cell_types <- prop[ ,rev (order (prop)) ]

cell_type1 <- paste (names (cell_types)[1] , cell_types [1], sep="=")
cell_type2 <- paste (names (cell_types)[2] , cell_types [2], sep="=")
cell_type3 <- paste (names (cell_types)[3] , cell_types [3], sep="=")
prop <- data.frame (cell_type1, cell_type2, cell_type3)
row.names (prop) <- barcodes[i]
propl[[i]] <- prop
}

prop <- do.call ("rbind", propl)

write.table (prop, "cell_proportions.txt", sep="\t", quote=F, row.names=TRUE, col.names=NA)



# Normalize per spot weights so cell type probabilities sum to 1 for each spot
# the normalized weights can now be considered as probabilities

weights <- list ()
for (i in (1:length (myRCTD@results))) {
x <- t (data.frame (myRCTD@results[[i]]$all_weights))
row.names (x) <- colnames(myRCTD@spatialRNA@counts) [i]
weights [[i]] <-x
}

weights <- do.call ("rbind", weights)
norm_weights <- normalize_weights (weights)
head (norm_weights)

cell_type_names <- colnames(norm_weights) # List of cell types

# Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
# Save each plot as a jpg file
for(i in 1:length(cell_type_names)){
   print (cell_type_names[i])
   plot_puck_continuous(myRCTD@spatialRNA, barcodes=colnames (myRCTD@spatialRNA@counts), norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=3)
   # ggsave(paste(resultsdir, cell_type_names[i],'_weights.jpg', sep=''), height=5, width=5, units='in', dpi=300)
}









