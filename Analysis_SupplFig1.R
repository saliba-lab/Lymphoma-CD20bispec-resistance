# R script for analysis of 15 biopsy samples of 7 lymphoma patients
# 
# related to Supplementary Figure 1.
# 
# seed for reproducibility
seed <- 1999

# load Seurat object.
ds <- readRDS("path/to/seurat_object")

# remove Erythrocyte genes HBB, HBA1, HBA2 from matrix.
# these genes were confounding factors due to RBC lysis issues
# in samples 6.1, 7.1 and 7.2
ds <- subset(ds, features = rownames(ds)[-c(which(rownames(ds) == to_ids["HBA1"]),
                                            which(rownames(ds) == to_ids["HBA2"]),
                                            which(rownames(ds) == to_ids["HBB"]))])

# perform quality control based on thresholds for
# UMI-counts, number of genes detected and fraction of mitochondiral counts

# calculate fraction mitochondrial counts
mt_counts <- Matrix::colSums(ds@assays$RNA@counts[
  features$ENSEMBL[which(stringr::str_detect(features$SYMBOL, "^MT-"))]
  , ]) 
ds@meta.data$percent.mt <- round(
  mt_counts / ds@meta.data$nCount_RNA * 100, 1
)

# get thresholds from Supplementary Table 3, Page 2
thresholds <- readxl::read_excel("path/to/SupplementaryTable_3.xlsx", sheet = "2")

# get cells with sufficient quality
cells <- c()
for (i in thresholds$Identifier) {
  cellsi <- colnames(subset(ds, subset = Identifier == i &
                              nCount_RNA > thresholds[i,2] &
                              nCount_RNA < thresholds[i,3] &
                              nFeature_RNA > thresholds[i,4] &
                              nFeature_RNA < thresholds[i,5] &
                              percent.mt > thresholds[i,6] &
                              percent.mt < thresholds[i,7]))

  cells <- c(cells, cellsi)
}

# filter out low_quality cells
ds <- subset(ds, cells = cells)

# Normalize mRNA counts, Identify variable features & scale
ds <- Seurat::NormalizeData(
  object               = ds,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
ds <- Seurat::FindVariableFeatures(
  object           = ds,
  selection.method = "vst", 
  nfeatures        = 3000
)

ds <- Seurat::ScaleData(ds)

# Dimensional reduction

ds <- Seurat::RunPCA(ds)

ds <- Seurat::RunUMAP(
  object = ds,
  dims = 1:50,
  seed.use = seed
)

# Clustering
set.seed(seed = seed)
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:50,
  reduction = "pca"
)

ds <- Seurat::FindClusters(
  object     = ds,
  resolution = 0.2
)

# Gene set scores for annotation
# dataframe of genesets (from Supplementary Table 3, Page 3)
genesets <- list(
  Bcell = to_ids[c("CD79A", "CD79B")],
  Tcell = to_ids[c("CD3D", "CD3G", "CD3E", "CD2", "IL32")],
  APC = to_ids[c("CST3", "LYZ", "TYROBP", "FCER1G", "APOC1", "SOD2")],
  Fibroblast = to_ids[c("COL3A1", "COL1A2", "MFAP5")],
  pDC = to_ids[c("CLEC4C", "IL3RA", "LILRB4", "SCT")]
)
# compute module scores
ds <- Seurat::AddModuleScore(
  object = ds, features = genesets, name = names(genesets), nbin = 24
)

# Cluster annotation based on module scores
# Cluster Annotations----------------------------------------------------------
Cluster.Annotation <- c(
  "0" = "B/malignant",
  "1" = "T cell",
  "2" = "B/malignant",
  "3" = "B/malignant",
  "4" = "B/malignant",
  "5" = "T cell",
  "6" = "B/malignant",
  "7" = "B/malignant",
  "8" = "B/malignant",
  "9" = "B/malignant",
  "10" = "T cell",
  "11" = "T cell",
  "12" = "T cell",
  "13" = "B/malignant",
  "14" = "B/malignant",
  "15" = "B/malignant",
  "16" = "B/malignant",
  "17" = "B/malignant",
  "18" = "Myeloid",
  "19" = "pDC",
  "20" = "Fibroblast"
)
ds@meta.data$Celltype <- factor(Cluster.Annotation[ds@meta.data$seurat_clusters],
                                unique(Cluster.Annotation)[c(1,2,3,4,5)])

# Calculation of cell type proportions
td2pct <- function(td) {
  for (i in 1:nrow(td)) {
    pcts <- c()
    for (ii in 1:ncol(td)) {
      pcts <- c(pcts, td[i,ii]/sum(td[i,]))
    }
    td[i, ] <- pcts*100
  }
  return(td)
}
proportions <- td2pct(table(ds@meta.data$Identifier, ds@meta.data$Celltype))
