## From url: https://github.com/satijalab/seurat/pull/10028
#devtools::install_github("satijalab/seurat-object", ref = "spaceranger-4.0", force = TRUE)
#devtools::install_github("satijalab/seurat", ref = "spaceranger-4.0", force = TRUE)

library(tidyverse)
library(Seurat)
library(reticulate)
library(anndata)
library(tidyquant)
library(nichenetr)
library(clusterProfiler)

#source(file.path("R","functions.R"))

python_dir="/Users/romansankowski/anaconda3/bin/python"
use_python(python_dir, required = T)
sys = import("sys")

## load data
file_path <- file.path("results","wt_proseg.h5ad")

adata <- anndata::read_h5ad(file_path)
adata$var_names_make_unique()

## transpose counts matrix
counts <- t(as.matrix(adata$X))

cell_metadata <- adata$obs
rownames(cell_metadata) <- cell_metadata$original_cell_id

seurat_obj <- CreateSeuratObject(counts = counts) %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  AddMetaData(cell_metadata) 

## extract variable genes
var_genes <- VariableFeatures(seurat_obj)

## set up query object
## add x y coordinates
coords <- adata$obsm$spatial %>%
  as.data.frame()
rownames(coords) <- colnames(seurat_obj)

seurat_obj@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coords
)

## run normalization
seurat_obj <- seurat_obj[, seurat_obj$volume < 2000] %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:20) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(dims = 1:20)

## reference was downloaded from cellxgene Skin of body - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - 10x
## url https://datasets.cellxgene.cziscience.com/1b84385a-0d2e-42ae-8d5a-db6706931060.h5ad
reference <- readRDS(file.path("data","Anderson.rds"))
reference2 <- readRDS(file.path("data","Anderson.rds"))

reference <- UpdateSeuratObject(reference)

## extract  data and rename genes
counts_ref <- reference[["RNA"]]$counts
genes <-  rownames(counts_ref)
genes <- bitr(genes, fromType = "ENSEMBL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

genes2 <- genes[!duplicated(genes$SYMBOL),]

counts_ref <- counts_ref[genes2$ENSEMBL,]
rownames(counts_ref) <- genes2$SYMBOL

## subset for var genes
var_genes2 <- var_genes[var_genes %in% rownames(counts_ref)]
counts_ref <- counts_ref[var_genes2,]

## extract metadata
meta_ref <- reference[[]]

## normalize
reference <- CreateSeuratObject(counts_ref) %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  AddMetaData(meta_ref)

## setup query object
query <- CreateSeuratObject(counts = seurat_obj[["RNA"]]$counts[var_genes2,]) %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

## run macrophage reference 
mac_reference <- readRDS(file.path("data","mac_filter_merge.RDS"))
DefaultAssay(mac_reference) <- "RNA"
meta_mac <- mac_reference[[]]

var_genes3 <- var_genes2[var_genes2 %in% rownames(mac_reference)]
mac_reference <- CreateSeuratObject(counts = mac_reference[["RNA"]]$counts[var_genes3,]) %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  AddMetaData(meta_mac)

## label transfer adjusted from url: https://satijalab.org/seurat/articles/integration_mapping.html
anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = reference$cell_type, dims = 1:30)
seurat_obj <- AddMetaData(seurat_obj, metadata = predictions)

summary(seurat_obj$prediction.score.max)
min(seurat_obj$prediction.score.max) ## its 0.2240658

## subset for high confidence celltype calls and renormalize
seurat_obj <- seurat_obj[,seurat_obj$prediction.score.max > .5] %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() 

## separately align macrophages to the macrophage reference

## extract macrophages
mac_query <- query[,colnames(seurat_obj)[seurat_obj$predicted.id == "professional antigen presenting cell"]] %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() 

## summarize clusters
mac_reference$new_clusters <- gsub("(C|H|I|R)","",mac_reference$new_clusters)

anchors <- FindTransferAnchors(reference = mac_reference, query = mac_query, dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = mac_reference$new_clusters, dims = 1:30)
mac_query <- AddMetaData(mac_query, metadata = predictions)

min(mac_query$prediction.score.max) ## ist 0.2956835
summary(mac_query$prediction.score.max) 

## add celltypes
celltype_pred <- seurat_obj$predicted.id
celltype_pred[colnames(mac_query)] <- paste0("Mac_",mac_query$predicted.id)

seurat_obj$cell_type_pred <- celltype_pred

saveRDS(seurat_obj, file.path("results","wt_skin_visiumHD_proseg_with_azimuth_label_transfer.rds"))

seurat_obj <- readRDS(file.path("results","wt_skin_visiumHD_proseg_with_azimuth_label_transfer.rds"))
write.csv(seurat_obj[[]][,c("original_cell_id","cell_type_pred")], file.path("results","azimuth_celltype_predictions.csv"))
