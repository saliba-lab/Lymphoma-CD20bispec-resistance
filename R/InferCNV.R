# Script for identification of widespread copy number variations
# from single cell data using inferCNV
# for MS4A1-positiva and MS4A1-negative B/malignant cells in Sample 1.4
# related to Figure 4a-d & Supplementary Figure 8

# load seurat object of Sample 1.4
ds <- readRDS("path/to/Seurat/Object.rds")

ds@meta.data$Icnv <- "malignant_CD20loss"
ds@meta.data$Icnv[ds@assays$RNA@counts[to_ids["MS4A1"],] > 1 & ds@meta.data$Celltype == "B/malignant"] <- "malignant_CD20pos"
ds@meta.data$Icnv[ds@meta.data$Celltype == "T cell"] <- "Normal_T"
ds@meta.data$Icnv[ds@meta.data$Celltype == "Myeloid"] <- "Normal_Myeloid"
write.table(data.frame(colnames(ds),
                       ds@meta.data$Icnv),
            "path/to/Annotation.txt", sep = "\t",
            col.names = F, row.names = F)
# make raw count matrix 
rcm <- Seurat::GetAssayData(ds, slot = "counts")

# create inferCNV object
infercnvObj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = rcm,
  annotations_file = "path/to/Annotation.txt",
  delim = "\t",
  gene_order_file = "path/to/gene_order_file.txt",
  ref_group_names = c("Normal_T", "Normal_Myeloid"),
  max_cells_per_group = 554
)

# run inferCNV
infercnvObj <- infercnv::run(
  infercnvObj,
  cutoff = 0.1,
  out_dir = "path/to/inferCNV/output/directory",
  cluster_by_groups = T,
  denoise = T,
  HMM = T
)
