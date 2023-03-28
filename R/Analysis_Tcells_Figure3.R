# Script for analysis of T cells of samples during CD20 bispec therapy
# excluding samples with less than 50 cells
# related to Figure 3 & Supplementary Figure 4
#
#

# load Seurat object containing T cells from samples during CD20 bispec therapy
# -> Samples 1.2, 1.3, 2.4, 2.5, 3.2, 3.3, 7.1, 7.2
ds <- readRDS("path/to/Seurat/object.rds")

# exclude samples contributing less than 50 T cells (Samples 1.2 & 1.3)
ds <- subset(ds, subset = Identifier %in% c("2.4", "2.5", "3.2", "3.3", "7.1", "7.2"))

# Randomly downsample remaining samples to a maximum of 1500 cells
# to keep the influence of individual samples onto the embedding and clustering in check


# Seurat analysis with bacth correction using fastMNN
ds <- Seurat::NormalizeData(
  object               = ds,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 10000
)
ds <- Seurat::FindVariableFeatures(
  object           = ds,
  selection.method = "vst", 
  nfeatures        = 1000
)

# MNN-corrected PCA, based on Normalized data 
set.seed(seed = seed)
ds@reductions[["mnn"]] <- Seurat::CreateDimReducObject(
  embeddings = batchelor::fastMNN(
    ds@assays$RNA@data[ds@assays$RNA@var.features, ],
    batch = ds@meta.data$sample,
    d     = 50
  )@int_colData$reducedDims$corrected,
  key        = "MNN_",
  assay      = "RNA"
)

# UMAP with MNN-corrected PCA
ds <- Seurat::RunUMAP(
  object         = ds,
  dims           = 1:50,
  seed.use       = seed,
  reduction      = "mnn",
  reduction.name = "mnn_umap",
  n.neighbors    = 50,
  min.dist       = 0.3,
  return.model = TRUE
  )

# Clustering
set.seed(seed = seed)
ds <- Seurat::FindNeighbors(
  object    = ds,
  dims      = 1:10,
  reduction = "mnn"
)

ds <- Seurat::FindClusters(
  object     = ds,
  resolution = 1
)

# Calculate module scores for gene sets to infer cluster identity
genesets <- list(
  Naive = to_ids[c("CCR7", "TCF7", "S1PR1", "LEF1")],
  Th = to_ids[c("CD4", "CD40LG")],
  Treg = to_ids[c("FOXP3", "IL2RA", "TNFRSF18", "TNFRSF4")],
  Th17 = to_ids[c("IL17A", "RORC", "CCR6", "CCL20", "IL4I1", "IL23R")],
  Tfh = to_ids[c("CD4", "CXCL13", "IL21")],
  CD8eff = to_ids[c("CD8A", "CD8B", "GZMA", "GZMB", "IFNG", "PRF1", "CCL4", "CCL5")],
  NKT = to_ids[c("CD8A", "NKG7", "FGFBP2", "ZNF683")],
  NK = to_ids[c("NKG7", "KLRF1", "TYROBP", "GNLY", "FCGR3A", "TRDC", "FCER1G", "NCAM1")],
  IFNresp = to_ids[c("ISG15", "RSAD2", "MX1", "IFIT3", "IFIT1", "OAS1")]
)

ds <- Seurat::AddModuleScore(
  object = ds, features = genesets, name = names(genesets), nbin = 24
)

# Annotate clusters
T.Annotation <- c(
  "0" = "CD8_eff",
  "1" = "CD4_Th",
  "2" = "CD8_naive",
  "3" = "CD8_eff",
  "4" = "CD8_eff",
  "5" = "CD4_Treg",
  "6" = "CD4_naive",
  "7" = "CD8_NKT",
  "8" = "CD4_Tfh",
  "9" = "NK",
  "10" = "IFN_responsive"
)
ds@meta.data$Annotation <- factor(T.Annotation[ds@meta.data$seurat_clusters],
                                  unique(T.Annotation)[c(1, 6, 8, 3, 5, 2, 7, 4, 9)])

# calculation of Exhaustion score for 3 CD8 cluster's cells
genesets <- list(Exhaustion = c("LAG3", "HAVCR2", "PDCD1", "PTMS", "FAM3C", "IFNG", "AKAP5", "CD7", "PHLDA1", "ENTPD1",
                        "SNAP47", "TNS3", "CXCL13","RDH10", "DGKH", "KIR2DL4", "LYST", "MIR155HG", "RAB27A",
                        "CSF1", "CTLA4", "TNFRSF9", "CD27", "CCL3", "ITGAE", "PAG1", "TNFRSF1B",
                        "GALNT1", "GBP2", "MYO7A", "ID3", "TOX", "ZBED2", "RBPJ", "ETV1", "NR4A2",
                        "ID2", "TRPS1", "PRDM1", "MAF", "EOMES", "VDR", "IKZF3", "MBD2", "BATF",
                        "EGR2", "BHLHE40", "TGIF1", "KMT2A"))

ds_CD8 <- Seurat::AddModuleScore(
  object = subset(ds, subset = Annotation2 %in% c("CD8_eff", "CD8_NKT", "CD8_naive")),
  features = genesets, name = names(genesets), nbin = 24, seed = 1999
)
