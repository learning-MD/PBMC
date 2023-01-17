library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)

pbmc <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR_undefined_1.12.23.rds")
pbmc@meta.data$cluster <- Idents(pbmc)
pbmc@meta.data$cluster <- as.character(pbmc@meta.data$cluster)

# Creating overall PBMC h5ad file for scCODA analysis in Python
SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat")
Convert("pbmc.h5Seurat", dest = "h5ad", assay = "RNA")

# CD4 T cell h5ad file
cd4 <- readRDS("PBMC_CellRanger_6.1.2/Upload/cd4_subclusters_1.12.23.rds")
cd4@meta.data$cluster <- Idents(cd4)
cd4@meta.data$cluster <- as.character(cd4@meta.data$cluster)

SaveH5Seurat(cd4, filename = "cd4.h5Seurat")
Convert("cd4.h5Seurat", dest = "h5ad", assay = "RNA")

# CD8 T cell h5ad file
cd8 <- readRDS("PBMC_CellRanger_6.1.2/Upload/cd8_subclusters_1.12.23.rds")
cd8@meta.data$cluster <- Idents(cd8)
cd8@meta.data$cluster <- as.character(cd8@meta.data$cluster)

SaveH5Seurat(cd8, filename = "cd8.h5Seurat")
Convert("cd8.h5Seurat", dest = "h5ad", assay = "RNA")

# NK cell h5ad file
nk <- readRDS("PBMC_CellRanger_6.1.2/Upload/NK_subclusters_1.12.23.rds")
nk@meta.data$cluster <- Idents(nk)
nk@meta.data$cluster <- as.character(nk@meta.data$cluster)

SaveH5Seurat(nk, filename = "nk.h5Seurat")
Convert("nk.h5Seurat", dest = "h5ad", assay = "RNA")
