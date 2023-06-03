library (openxlsx)
library (spacexr)
library (Seurat)
library (dplyr)
library (ggplot2)

myRCTD <- readRDS ("myRCTD_visium_allen_full_annot_multi_mode.rds")

weights <- list ()
for (i in (1:length (myRCTD@results))) {
x <- t (data.frame (myRCTD@results[[i]]$all_weights))
row.names (x) <- colnames(myRCTD@spatialRNA@counts) [i]
weights [[i]] <-x
}

weights <- do.call ("rbind", weights)
norm_weights <- normalize_weights (weights)
head (norm_weights)


## Add the manual annotation

annot <- read.delim ("/Volumes/texas/iit_projects/martina/Northwestern University/NUSeq Core Facility - Martina03_9.16.2021/WORKING/Location information/G1_1A_shmv01.txt")
colnames (annot) [1] <- "Barcode"
colnames (annot) [3] <- "Selection"
head (annot)

idx <- match (row.names (norm_weights), annot$Barcode)
norm_weights <- cbind (norm_weights, data.frame (selection= annot$Selection[idx]))
norm_weights$selection[is.na (norm_weights$selection)] <- 0
norm_weights$selection[norm_weights$selection == ""] <- 0
table (norm_weights$selection)

norm_weights$selection[norm_weights$selection == "CA1"] <- 1
norm_weights$selection[norm_weights$selection == "CA2"] <- 0.2
norm_weights$selection[norm_weights$selection == "CA3"] <- 0.4
norm_weights$selection[norm_weights$selection == "DG"] <- 0.6

norm_weights$selection <- as.numeric (norm_weights$selection)
table (norm_weights$selection)

cell_type_names <- colnames(norm_weights) # List of cell types


vals <- as.numeric (norm_weights$selection)
names (vals) <- row.names (norm_weights)

# Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
# Save each plot as a jpg file
midpoint <- 0.4
i <- 42
   print (cell_type_names[i])
   plot_puck_continuous(myRCTD@spatialRNA, barcodes=colnames (myRCTD@spatialRNA@counts), vals, title =cell_type_names[i], size=0.3) + ggplot2::scale_colour_gradient2(midpoint = midpoint, low="white", mid="blue", high="darkred")
   ggsave('preselected_cells.jpg', height=8, width=8, units='in', dpi=300)









