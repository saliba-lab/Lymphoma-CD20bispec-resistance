# Script for analysis of Samples 3.2 & 3.3
# related to Figure 4f and Supplementary Figure 10

# load Seurat object including cells from Samples 3.2 & 3.3
ds <- readRDS("path/to/Seurat/Object.rds")

# basic Seurat workflow
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 1000)
ds <- Seurat::ScaleData(ds)
ds <- Seurat::RunPCA(ds)
ds <- Seurat::RunUMAP(ds, dims = 1:10)
ds <- Seurat::FindNeighbors(ds, dims = 1:10)
ds <- Seurat::FindClusters(ds, resolution = 0.1)
