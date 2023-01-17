library(Seurat)
library(sctransform)
library(harmony)
library(patchwork)
library(tidyverse)
library(Nebulosa)

# Loading data
pbmc.combined <- readRDS("PBMC_CellRanger_6.1.2/Upload/decontX_all_combined.rds")

# Including TCR/BCR genes
# SCTransform individual batches (since 5 samples were multiplexed in each batch)
features <- SplitObject(pbmc.combined, split.by = "batch")

features <- lapply(X = features,
                   FUN = SCTransform,
                   vst.flavor = "v2",
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)

pbmc <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)

VariableFeatures(pbmc) <- var.features

saveRDS(pbmc, "PBMC_CellRanger_6.1.2/Upload/pbmc_including_TCR_BCR_postSCT.rds")
rm(list = ls())
gc()

pbmc <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_including_TCR_BCR_postSCT.rds")

# Batch correction with Harmony
pbmc_correct <- RunPCA(pbmc, verbose = TRUE)
pbmc_correct <- RunHarmony(pbmc_correct, assay.use="SCT", group.by.vars = "batch")
pbmc_correct <- RunUMAP(pbmc_correct, reduction = "harmony", dims = 1:30)
pbmc_correct <- FindNeighbors(pbmc_correct, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.4)

DimPlot(pbmc_correct, label = TRUE)
saveRDS(pbmc_correct, "PBMC_CellRanger_6.1.2/Upload/pbmc_post_processing_include_TCR_BCR.rds")
ggsave("PBMC_CellRanger_6.1.2/Upload/Include_TCR_dimplot.png", dpi = 300)

# Excluding TCR/BCR genes
rm(list = ls())
gc()

pbmc.combined <- readRDS("PBMC_CellRanger_6.1.2/Upload/decontX_all_combined.rds")

# Removing TCR genes
pbmc.combined <- pbmc.combined[!grepl("^TR[ABDG][VJC]", rownames(pbmc.combined)), ]

# Removing BCR-genes
pbmc.combined <- pbmc.combined[!grepl("^IG[HKL]V", rownames(pbmc.combined)), ]
pbmc.combined <- pbmc.combined[!grepl("^IG[HKL]J", rownames(pbmc.combined)), ]
pbmc.combined <- pbmc.combined[!grepl("^IG[KL]C", rownames(pbmc.combined)), ]
pbmc.combined <- pbmc.combined[!grepl("^IGH[ADEGM]", rownames(pbmc.combined)), ]

# SCTransform individual batches
features <- SplitObject(pbmc.combined, split.by = "batch")

features <- lapply(X = features,
                   FUN = SCTransform,
                   vst.flavor = "v2",
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)

pbmc <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)

VariableFeatures(pbmc) <- var.features

saveRDS(pbmc, "PBMC_CellRanger_6.1.2/Upload/pbmc_excluding_TCR_BCR_postSCT.rds")
rm(list = ls())
gc()

pbmc <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_excluding_TCR_BCR_postSCT.rds")

# Batch correction with Harmony
pbmc_correct <- RunPCA(pbmc, verbose = TRUE)
pbmc_correct <- RunHarmony(pbmc_correct, assay.use="SCT", group.by.vars = "batch")
pbmc_correct <- RunUMAP(pbmc_correct, reduction = "harmony", dims = 1:30)
pbmc_correct <- FindNeighbors(pbmc_correct, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.4)

DimPlot(pbmc_correct, label = TRUE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Exclude_TCR_dimplot.png", dpi = 300)
saveRDS(pbmc_correct, "PBMC_CellRanger_6.1.2/Upload/pbmc_post_processing_exclude_TCR_BCR.rds")

pbmc_correct <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_post_processing_exclude_TCR_BCR.rds")

# Finding cluster markers
cluster.markers.wilcox_only.pos <- FindAllMarkers(pbmc_correct, assay = "RNA",
                                                  logfc.threshold = 0.25,
                                                  min.pct = 0.1,
                                                  min.diff.pct = 0.25,
                                                  test.use = "wilcox",
                                                  only.pos = TRUE)

write_csv(cluster.markers.wilcox_only.pos, "PBMC_CellRanger_6.1.2/Upload/cluster_markers_wilcox.csv")

# Visualize distribution of markers

# Broad T-cell markers
FeaturePlot(pbmc_correct, features = c("CD3E", "CD3D", "CD8A", "CD8B"))
plot_density(pbmc_correct, features = c("CD3E", "CD3D", "CD8A", "CD8B"))

# Cytotoxic T cells
FeaturePlot(pbmc_correct, features = c("CD3D", "GZMK", "GZMB"))
plot_density(pbmc_correct, features = c("CD3D", "GZMK", "GZMB"))

# Tregs
FeaturePlot(pbmc_correct, features = c("CD3D", "FOXP3", "CTLA4"))
plot_density(pbmc_correct, features = c("CD3D", "FOXP3", "CTLA4"))

# NK cell markers
FeaturePlot(pbmc_correct, features = c("GNLY", "NKG7", "KLRD1"))
plot_density(pbmc_correct, features = c("GNLY", "NKG7"))

# B cell markers
FeaturePlot(pbmc_correct, features = c("MS4A1", "CD79A"))

# Monocyte markers
FeaturePlot(pbmc_correct, features = c("CD14", "FCGR3A"))

# Dendritic markers
FeaturePlot(pbmc_correct, features = c("FCER1A", "CD1C", "FSCN1",
                                       "ITGAX", "CD14", "CLEC7A",
                                       "CLEC6A", "IL3RA", "CX3CR1",
                                       "LILRA4", "CST3", "THBD"))

FeaturePlot(pbmc_correct, features = c("CLEC9A", "IRF8",
                                       "CLEC10A", "FCER1A", "CD1C",
                                       "FSCN1", "CCR7"))

FeaturePlot(pbmc_correct, features = c("HLA-DQA1", "HLA-DQB1", "MS4A7"))

# Platelets
FeaturePlot(pbmc_correct, features = c("PPBP"))

# Annotation
cluster.ids <- c("CD4+ T cells",
                 "CD14+ monocytes I",
                 "CD14+ monocytes II",
                 "CD8+ T cells I",
                 "NK cells I",
                 "CD8+ T cells II",
                 "CD16+ monocytes",
                 "B cells",
                 "CD8+ T cells III",
                 "NK cells II",
                 "Plasmacytoid dendritic cells I",
                 "CD1c+ dendritic cells",
                 "Platelets",
                 "Plasmacytoid dendritic cells II",
                 "Undefined"
)

names(cluster.ids) <- levels(pbmc_correct)
pbmc_correct <- RenameIdents(pbmc_correct, cluster.ids)

DimPlot(pbmc_correct, label = TRUE) + NoLegend() + ggtitle("Annotated Clusters")
ggsave("PBMC_CellRanger_6.1.2/Upload/Exclude_TCR_DimPlot_annotated.png", dpi = 300)

saveRDS(pbmc_correct, "PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR.rds")
