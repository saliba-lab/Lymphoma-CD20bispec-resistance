# Script for analysis of Sample 5.1 T cells
# related to Supplementary Figure 6

# load Seurat object with T cells from Sample 5.1
ds <- readRDS("path/to/Seurat/object.rds")

# basic Seurat workflow
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 1000)
ds <- Seurat::ScaleData(ds)
ds <- Seurat::RunPCA(ds)
ds <- Seurat::RunUMAP(ds, dims = 1:10)
ds <- Seurat::FindNeighbors(ds, dims = 1:10)
ds <- Seurat::FindClusters(ds, resolution = 2.5)

# Annotation of T cell subpopulations
T.Annotation <- c(
  "0" = "CD4 Th",
  "1" = "CD8 eff",
  "2" = "Prolif",
  "3" = "CD4 Th",
  "4" = "CD4 naive",
  "5" = "CD8 eff",
  "6" = "CD4 naive",
  "7" = "CD8 eff",
  "8" = "CD4 Th",
  "9" = "CD4 Th",
  "10" = "CD4 naive",
  "11" = "NK",
  "12" = "NK",
  "13" = "Treg",
  "14" = "CD8 eff"
)
