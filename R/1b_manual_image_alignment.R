library(STutility)
library(Seurat)
library(tidyverse)

seurat_obj <- readRDS(file.path("results","wt_skin_visiumHD_proseg_with_azimuth_label_transfer.rds"))

ManualAlignImages(
  seurat_obj,
  type = "masked.masks",
  reference.index = 1,
  edges = TRUE,
  verbose = FALSE,
  limit = 0.3,
  maxnum = 1000,
  fix.axes = FALSE,
  custom.edge.detector = NULL
)
