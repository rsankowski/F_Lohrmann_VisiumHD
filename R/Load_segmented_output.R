library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)


install.packages("SeuratObject")
devtools::install_github("satijalab/seurat-object", ref = "spaceranger-4.0", force = TRUE)
devtools::install_github("satijalab/seurat", ref = "spaceranger-4.0", force = TRUE)

library(SeuratObject)
library(Seurat)

## load segmented output
data.dir <- file.path("/Users/romansankowski/Documents/single_cell_analysis/F_Lohrmann_VisiumHD/data/WT2_spaceranger4/outs")
v3d_polys_only <- Load10X_Spatial(data.dir = data.dir,
                                  slice = "slice1",
                                  bin.size = c("polygons"))
