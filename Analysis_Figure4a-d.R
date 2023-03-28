# Script for Analysis of CD20 recovery in patient 1 sample 1.4
# related to Figure 4a-d and Supplementary Figure 8

# load Seurat object containing Sample 1.4 B/malignant cells
ds <- readRDS("path/to/Seurat/Object/Sample1.4_Bcells.Rds")

# basic Seurat analysis
ds <- Seurat::NormalizeData(ds)
ds <- Seurat::FindVariableFeatures(ds, nfeatures = 1000)
ds <- Seurat::ScaleData(ds)
ds <- Seurat::RunPCA(ds)
ds <- Seurat::FindNeighbors(ds, dims = 1:10)
ds <- Seurat::FindClusters(ds, resolution = 0.05)
ds <- Seurat::RunUMAP(ds, dims = 1:10)

# make meta.data column with information if MS4A1 is expressed
# which is the case when MS4A1 expression is greater than 0
ds@meta.data$CD20 <- "neg"
ds@meta.data$CD20[ds@assays$RNA@data[to_ids["MS4A1"],] > 0] <- "pos"

# DE gene expression between CD20 positive and CD20 negative -> not informative
ds <- Seurat::SetIdent(ds, value = ds@meta.data$CD20)
markers <- Seurat::FindMarkers(ds, ident.1 = "pos", ident.2 = "neg", only.pos = F)

# clusterProfiler GO enrichment
# select upregulated genes with p-value less than 0.05
CP_markers <- subset(markers, subset = markers$avg_log2FC > 0 & markers$p_val < 0.05)
CP_up <- rownames(CP_markers)
# GO enrichment
CP_GO <- clusterProfiler::enrichGO(CP_up, # vector of up regulated genes
                                       "org.Hs.eg.db", # orgdb= package that contains gene label types correspondances
                                       keyType = "ENSEMBL", # indicate that genes are labeled using symbols
                                       ont = "BP")
# use simplify to get rid of redundant GO terms
CP_GO <- clusterProfiler::simplify(CP_GO)

# Identify actively cycling cells
# get cell cycle signature genes for s-phase and g2m-phase
cc_file <- RCurl::getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv")
cell_cycle_genes <- read.csv(text = cc_file)
s.genes <- cell_cycle_genes$geneID[cell_cycle_genes$phase == "S"]
g2m.genes <- cell_cycle_genes$geneID[cell_cycle_genes$phase == "G2/M"]

# Score cell cycle signatures across cells
ds <- Seurat::CellCycleScoring(ds,
                                  s.features = s.genes,
                                  g2m.features = g2m.genes)
