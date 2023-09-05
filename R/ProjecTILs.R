# Script for Projecting T-cells with ProjecTILs

options(timeout = max(900, getOption("timeout")))

download.file("https://figshare.com/ndownloader/files/38921366", destfile = "/Path/to/CD8T_human_ref_v1.rds")
ref.cd8 <- ProjecTILs::load.reference.map("/Path/to/CD8T_human_ref_v1.rds")

ds <- readRDS("/Path/to/CD8_effector_T_cell_object.Rds")

ref.cd8@meta.data$species <- "human"
ds@meta.data$species <- "human"

ds <- ProjecTILs::Run.ProjecTILs(ds, ref = ref.cd8, filter.cells = FALSE)

ProjecTILs::plot.projection(ref.cd8, ds,
                            linesize = 0.3, pointsize = 1, ref.alpha = 0.9, ref.size =  1)
