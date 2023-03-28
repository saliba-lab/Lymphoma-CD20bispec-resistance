# SoupX ambient RNA removal

filedir <- "path/to/datasets/"

datasets <- c("G3", "G4", "G5", "H4", "H5",
              "H6", "H11", "A1_2", "A9_2", "C2_2",
              "C3_2", "C4_2", "C10_2", "E7_3", "E8_3")

for (dataset in datasets){
  # Create SoupChannel object
  print(dataset)
  from_filt <- paste0(filedir, dataset, "/filtered_feature_bc_matrix/")
  from_raw <- paste0(filedir, dataset, "/raw_feature_bc_matrix/")
  toc <- Seurat::Read10X(from_filt)  
  tod <- Seurat::Read10X(from_raw)  
  sc <- SoupX::SoupChannel(tod, toc) 
  
  # Create Seurat object and make clustering for SoupX
  sd <- Seurat::CreateSeuratObject(counts = toc)
  sd <- Seurat::NormalizeData(sd, verbose = F)
  sd <- Seurat::FindVariableFeatures(sd)
  sd <- Seurat::ScaleData(sd)
  sd <- Seurat::RunPCA(sd, verbose = F)
  sd <- Seurat::RunUMAP(sd, dims = 1:30, verbose = F)
  sd <- Seurat::FindNeighbors(sd, dims = 1:30, verbose = F)
  sd <- Seurat::FindClusters(sd, resolution = 0.3, verbose = T) 
  
  # Set the SoupChannel clusters
  sc <- SoupX::setClusters(sc, sd@meta.data$seurat_clusters)
  
  # ambient RNA estimation by SoupX. With very homogeneous data this will not work
  # therefore estimations above 0.3 will be regulated to 0.05 (middle contamination)
  sc <- SoupX::autoEstCont(sc, forceAccept = TRUE)
  print(mean(sc$metaData$rho))
  if (mean(sc$metaData$rho) > 0.3) {
    sc <- SoupX::setContaminationFraction(sc, 0.05)
  }
  
  # adjust matrix 
  adj.matrix  <- SoupX::adjustCounts(sc, roundToInt = T)
  to <- paste0(filedir, dataset, "/SoupX_matrix.mtx")
  Matrix::writeMM(obj = adj.matrix, file = to)
}
