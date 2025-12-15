library(tidyverse)
library(ggthemes)
library(viridis)
library(patchwork)
library(tessera)

## set parameters
fig.size <- function(h, w) {
  options(repr.plot.height = h, repr.plot.width = w)
}

verbose = TRUE
show_plots = TRUE

###### STEP 0 ######
npcs = 10
## Graph pruning
prune_thresh_quantile = 0.95
prune_min_cells = 5


###### STEP 1: GRADIENTS ######
smooth_distance = c('none', 'euclidean', 'projected', 'constant')[3] 
smooth_similarity = c('none', 'euclidean', 'projected', 'constant')[3] 


###### STEP 2: DMT ######
## ... no options


###### STEP 3: AGGREGATION ######
max_npts = 50
min_npts = 5
alpha = 1 ## 0.2 = conservative merging, 2 = liberal merging 

## load data

# Read AnnData file
file_path <- file.path("results","wt_skin_visiumHD_proseg_with_azimuth_label_transfer.rds")

seurat_obj <- readRDS(file_path)

counts <- seurat_obj[["RNA"]]$counts
meta_data <- GetTissueCoordinates(seurat_obj)[,-c(3:5)]
colnames(meta_data)[1:2] <- c("X","Y")

##add celltype info
meta_data$type <- seurat_obj$cell_type_pred

## transpose counts
counts <- t(as.matrix(counts)) 

## plot metadata
fig.size(8, 8)
ggplot() + 
  geom_point(data = meta_data, aes(X, Y, color = type), size=.5) + 
  theme_void() + 
  scale_color_tableau(palette = "Tableau 20") + 
  coord_sf(expand = FALSE) + 
  NULL

## Prepare data
meta_vars_include = c('type')
dmt = init_data(meta_data$X,meta_data$Y, counts, meta_data, meta_vars_include)
dmt = prune_graph(dmt, thresh_quantile = prune_thresh_quantile, mincells = prune_min_cells) 

## add triangles
dmt = add_exterior_triangles(dmt)

fig.size(80, 80)

if (show_plots) {    
  ggplot() + 
  
    ## small data 
    geom_segment(data = dmt$edges, aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'black', lwd = .2) + 
    geom_segment(data = dmt$edges[boundary == TRUE], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'red', lwd = .2) + 
    geom_point(data = dmt$pts, aes(X, Y), size = .5) + 
    
    theme_void(base_size = 20) + 
    coord_cartesian(expand = FALSE) + 
    labs(title = 'Pruned adjacency graph') + 
    NULL
}

## pca
dmt$udv_cells = do_pca(dmt$counts, npcs)

# Step 1: compute gradients on all data structures
field = compute_gradients(dmt, smooth_distance, smooth_similarity)
field = compress_gradients_svd(field)

if (show_plots) {    
  
  len_plot_constant = .8
  fig.size(15, 15)
  # fig.size(20, 20)
  ggplot() + 
    geom_point(data = dmt$pts, aes(X, Y), size = .5) + 
    geom_segment(data = dmt$edges[boundary == TRUE, ], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'red') + 
    
    ## Triangle Gradients
    geom_segment(
      data = data.table(dmt$tris, field$tris_svd), 
      aes(
        x=X-len_plot_constant*(len_grad+len_ortho)*dx_ortho, 
        y=Y-len_plot_constant*(len_grad+len_ortho)*dy_ortho, 
        xend=X+len_plot_constant*(len_grad+len_ortho)*dx_ortho, 
        yend=Y+len_plot_constant*(len_grad+len_ortho)*dy_ortho
      ), 
      linewidth = .4, alpha = 1, 
      color = 'blue'
    ) + 
    
    theme_void() + 
    coord_fixed(expand = FALSE) + 
    NULL
}

# Step 2: DMT

## compute f
dmt = dmt_set_f(dmt, field)

if (show_plots) {    
  ntri = max(which(dmt$tris$external == FALSE))
  i = Matrix::t(dmt$tri_to_pt[1:ntri, ])@i+1
  plt_df = data.table(
    X = dmt$pts$X[i],
    Y = dmt$pts$Y[i],
    f = rep(dmt$tris$f[1:ntri], each = 3)
  )[
    , id := rep(1:ntri, each = 3)
  ][]
  
  
  fig.size(15, 15)
  ggplot() + 
    geom_polygon(data = plt_df, aes(X, Y, group = id, fill = f, color = f)) + 
    theme_void() + 
    coord_fixed(expand = FALSE) + 
    scale_fill_viridis() + 
    scale_color_viridis() + 
    NULL
}

## forests
dmt$prim = do_primary_forest(dmt)
dmt$dual = do_dual_forest(dmt)

if (show_plots) {    
  fig.size(15, 15)
  ggplot() +     
    ## primary forest
    geom_point(data = dmt$tris[dmt$dual$maxima, ], aes(X, Y), color = 'blue', size = 2) + 
    geom_segment(data = dmt$dual$edges, aes(x=x0, y=y0, xend=x1, yend=y1), color = 'blue') + 
    
    ## primary forest
    geom_point(data = dmt$pts[dmt$prim$minima, ], aes(X, Y), color = 'red', size = 2) + 
    geom_segment(data = dmt$prim$edges, aes(x=x0, y=y0, xend=x1, yend=y1), color = 'red') + 
    
    theme_void() + 
    coord_cartesian(expand = FALSE) + 
    NULL
}


## extract epaths
dmt$e_sep = dmt_get_separatrices(dmt)

if (show_plots) {    
  fig.size(15, 15)
  ggplot() + 
    
    geom_segment(data = dmt$edges[dmt$e_sep, ], aes(x = x0_tri, y = y0_tri, xend = x1_tri, yend = y1_tri), lwd = 1, color = 'blue') + 
    geom_segment(data = dmt$edges[boundary == TRUE], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'blue', lwd = 1) + 
    
    ## primary forest
    geom_point(data = dmt$pts[dmt$prim$minima, ], aes(X, Y), color = 'red', size = 2) + 
    geom_segment(data = dmt$prim$edges, aes(x=x0, y=y0, xend=x1, yend=y1), color = 'red') + 
    
    theme_void() + 
    coord_cartesian(expand = FALSE) + 
    NULL
}

## extract tiles
dmt = dmt_assign_tiles(dmt)
aggs = dmt_init_tiles(dmt)

if (show_plots) {    
  set.seed(2)
  fig.size(15, 20)
  ggplot() + 
    geom_sf(data = aggs$meta_data$shape) + 
    # geom_point(data = dmt$pts, aes(X, Y, color = factor(agg_id, sample(nrow(aggs$meta_data)))), size = 1) + 
    # scale_color_tableau() + 
    theme_void() + 
    coord_sf(expand = FALSE) + 
    # coord_cartesian(expand = FALSE) + 
    guides(color = 'none') + 
    
    NULL 
}

if (show_plots) {    
  set.seed(2)
  fig.size(15, 20)
  ggplot() + 
    geom_sf(data = aggs$meta_data$shape) + 
    geom_point(data = dmt$pts, aes(X, Y, color = type), size=.5) + 
    theme_void() + 
    coord_sf(expand = FALSE) + 
    scale_color_tableau() + 
    guides(color = 'none') + 
    NULL 
}

# Step 3: Aggregation

## Merge main
aggs = init_scores(aggs, agg_mode=2, alpha=alpha, max_npts=max_npts)
aggs = merge_aggs(aggs, agg_mode=2, max_npts=max_npts)
dmt = update_dmt_aggid(dmt, aggs)
aggs = update_agg_shapes(dmt, aggs)

## Merge small outliers
aggs = init_scores(aggs, agg_mode=3, alpha=alpha, min_npts=min_npts)
aggs = merge_aggs(aggs, agg_mode=3, min_npts=min_npts)
dmt = update_dmt_aggid(dmt, aggs)
aggs = update_agg_shapes(dmt, aggs)

## Final tiles 
if (show_plots) { 
  purrr::map(1:3, function(i) {
    ggplot(cbind(aggs$meta_data, val=aggs$pcs[, i])) + 
      geom_sf(aes(geometry = shape, fill = val)) + 
      theme_void(base_size = 16) + 
      coord_sf(expand = FALSE) + 
      scale_fill_gradient2_tableau() + 
      guides(color = 'none') + 
      labs(title = paste0('PC', i)) + 
      NULL 
  }) %>% 
    purrr::reduce(`|`)
}

# Results

## Aggregates 
## The primary output is the tiles. Each tile has a row in the meta_data table: 
##  - npts denotes the number of cells in the tile. 

head(aggs$meta_data)


## We also have pooled gene counts, for differential gene expression analysis. 
aggs$counts[1:5, 1:5]

## And we have PCA embeddings for the tiles. 
head(aggs$pcs)

## The rest of the fields are internal to the algorithm and can be ignored. 

setdiff(names(aggs), c('pcs', 'meta_data', 'counts'))

## Get tiles
res = GetTiles(
  X = meta_data$X, 
  Y = meta_data$Y, 
  counts = counts, 
  meta_data = meta_data, 
  meta_vars_include = meta_vars_include,
)
dmt = res$dmt
aggs = res$aggs

## plot
fig.size(10, 30)
purrr::map(1:3, function(i) {
  ggplot(cbind(aggs$meta_data, val=aggs$pcs[, i])) + 
    geom_sf(aes(geometry = shape, fill = val)) + 
    theme_void(base_size = 16) + 
    coord_sf(expand = FALSE) + 
    scale_fill_gradient2_tableau() + 
    guides(color = 'none') + 
    labs(title = paste0('PC', i)) + 
    NULL 
}) %>% 
  purrr::reduce(`|`)

## Cluster and label tiles 
obj = Seurat::CreateSeuratObject(
  counts = aggs$counts, 
  meta.data = tibble::column_to_rownames(data.frame(dplyr::select(aggs$meta_data, -shape)), 'id')
)

## Seurat doesn't do sf shapes well 
obj@meta.data$shape = aggs$meta_data$shape

## Represent each tile as the mean PC embeddings of all its cells 
## NOTE: this tends to produce more biologically meaningful results than pooling gene counts per tile 
rownames(aggs$pcs) = colnames(obj)
obj[['pca']] = Seurat::CreateDimReducObject(embeddings = aggs$pcs, loadings = dmt$udv_cells$loadings, key = 'pca_', assay = Seurat::DefaultAssay(obj))

## clustering
.verbose = FALSE
obj = obj %>% 
  NormalizeData(normalization.method = 'LogNormalize', scale.factor = median(obj@meta.data$nCount_RNA), verbose = .verbose) %>% 
  RunUMAP(verbose = .verbose, dims = 1:10, reduction = 'pca') %>% 
  Seurat::FindNeighbors(features = 1:10, reduction = 'pca', verbose = .verbose) %>% 
  Seurat::FindClusters(verbose = .verbose, resolution = c(2))

## Let's see the aggregate clusters in UMAP and physical space. 
p1 = DimPlot(obj, reduction = 'umap', group.by = 'seurat_clusters') + scale_color_tableau('Classic 10') 
p2 = ggplot(obj@meta.data) + 
  geom_sf(aes(geometry = shape, fill = seurat_clusters)) + 
  theme_void(base_size = 16) + 
  coord_sf(expand = FALSE) + 
  scale_fill_tableau('Classic 10') + 
  NULL 

fig.size(6, 12)
(p1 | p2) + plot_layout(widths = c(1, 1))

### Transfer agg information to cells
dmt$pts$spatial_cluster = obj@meta.data$seurat_clusters[dmt$pts$agg_id]

p1 = ggplot() + 
  geom_sf(data = obj@meta.data, aes(geometry = shape), fill = NA) + 
  geom_point(data = dmt$pts, aes(X, Y, color = type), size=.1) + 
  scale_color_tableau(palette = 'Tableau 20') + 
  theme_void() + 
  coord_sf(expand = FALSE) + 
  NULL
p2 = ggplot() + 
  geom_sf(data = obj@meta.data, aes(geometry = shape, fill = seurat_clusters), alpha = .2) + 
  geom_point(data = dmt$pts, aes(X, Y, color = spatial_cluster), size=.1) + 
  scale_color_tableau('Classic 10') + 
  scale_fill_tableau('Classic 10') + 
  theme_void() + 
  guides(fill = 'none') + 
  coord_sf(expand = FALSE) + 
  NULL
fig.size(16, 40)
p1 | p2

ggsave(file.path("plots","spatial","cell_types_tessera_clusters.pdf"), height = 10, width = 12)

## Let's look at the composition of the spatial clusters. 
fig.size(8, 10)
dmt$pts %>% 
    with(table(type, spatial_cluster)) %>% 
    prop.table(2) %>% 
    data.table() %>% 
    ggplot(aes(spatial_cluster, 100 * N, fill = type)) + 
        geom_bar(stat = 'identity', position = position_stack()) + 
        scale_fill_tableau(palette = 'Tableau 20') + 
        theme_bw(base_size = 20) + 
        labs(y = '% of spatial cluster', fill = 'cell type') + 
        NULL