## See https://www.10xgenomics.com/resources/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
## See https://github.com/dmcable/spacexr

#options(timeout = 600000000)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)


library (openxlsx)
library (spacexr)

# load the visium data
dir <- "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/Space Ranger Output/G1-1A"
visium <- read.VisiumSpatialRNA (dir)

# Create a list of barcodes from the column names of the count matrix
barcodes <- colnames(visium@counts)

# Plot number of UMIs per barcode (spot)
## The number of UMIs (RNA species) is proportional to the number of cells 

plot_puck_continuous(puck=visium, barcodes=barcodes, plot_val=visium@nUMI, 
                     size=1, ylimit=c(0,round(quantile(visium@nUMI,0.9))), 
                     title='plot of nUMI') 
           
           
           
## get an independent chromium single cell RNA-Seq data with some hippocampus annotated cell types
## this is the reference annotated single cell data

refdir <- system.file("extdata",'Reference/Visium_Ref',package = 'spacexr') # directory for the reference

counts <- read.csv(file.path(refdir,"counts.csv")) # load in cell types
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
dim (counts)
# 307 genes * 510 cell types

cell_types <- read.csv(file.path(refdir,"cell_types.csv")) # load in cell types
cell_types <- setNames(cell_types[[2]], cell_types[[1]])
cell_types <- as.factor(cell_types) # convert to factor data type
table (cell_types)

nUMI <- read.csv(file.path(refdir,"nUMI.csv")) # load in cell types
nUMI <- setNames(nUMI[[2]], nUMI[[1]])

reference <- Reference(counts, cell_types, nUMI)



## RCTD in multimode

Sys.setenv("OPENBLAS_NUM_THREADS"=2)

myRCTD <- create.RCTD(visium, reference, max_cores = 2)

# full mode, which assigns any number of cell types per spot and is recommended for technologies with poor spatial resolution such as 100-micron resolution Visium; 
# multi mode, an extension of doublet mode that can discover more than two cell types per spot as an alternative option to full mode

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
saveRDS (myRCTD, "myRCTD_visium_multi_mode.rds")



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


## Comparison with the manual annotation

annot <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01.txt")
colnames (annot) [1] <- "Barcode"
head (annot)

annot <- merge (annot, prop, by.x="Barcode", by.y= "row.names", all.x=TRUE)
write.xlsx (annot, "/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01_annot.xlsx")



