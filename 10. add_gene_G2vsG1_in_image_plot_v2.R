library (Seurat)         # Seurat_4.9.9.9041
library (ggpubr)
library (sctransform)    # sctransform_0.3.5

### Slide 1
brain <- readRDS ("brain_slide1_G1G3_groups.rds")

### Slide 2
brain2 <- readRDS ("brain_slide2_G2G4_groups.rds")

brain <- merge(brain, y = c(brain2))

DefaultAssay(brain) <- "Spatial"
brain <- JoinLayers (brain)

# Here we normalize everything with SCT
brain <- SCTransform(brain, vst.flavor = "v2") 

Images (brain)

## FIXME Continue
