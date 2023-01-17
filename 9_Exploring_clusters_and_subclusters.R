library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)
library(Nebulosa)
library(harmony)
library(SingleR)

# Load original PBMC dataset
pbmc <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR.rds")

# Remove the undefined group
pbmc <- subset(pbmc, idents = "Undefined", invert= TRUE)


# Defining the CD14+ monocytes II cluster as CD141+ dendritic cells
markers <- FindAllMarkers(pbmc,
                          assay = "RNA",
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          min.diff.pct = 0.25,
                          test.use = "roc",
                          only.pos = TRUE)

plot_density(pbmc, features = c("THBD", "CLEC9A"))

plot_density(pbmc, features = c("AXL", "SIGLEC6"))

cluster.ids <- c("CD4+ T cells",
                 "CD14+ monocytes",
                 "CD141+ dendritic cells",
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
                 "Plasmacytoid dendritic cells II"
)

names(cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cluster.ids)

DimPlot(pbmc, label = TRUE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/DimPlot_annotated_1.13.23.png", dpi = 600)

saveRDS(pbmc, "PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR_undefined_1.12.23.rds")

# Save plots
DimPlot(pbmc, label = FALSE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/DimPlot_legend.png", dpi = 600)

DimPlot(pbmc, label = TRUE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/DimPlot_legend_and_labels.png", dpi = 600)

DimPlot(pbmc, label = TRUE) + NoLegend()
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/DimPlot_label.png", dpi = 600)

# Save RDS file
saveRDS(pbmc, "PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR_undefined_1.12.23.rds")
rm(list = ls())
gc()
pbmc <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR_undefined_1.12.23.rds")

# Subcluster CD4+ T cells
cd4_t_cell <- subset(pbmc, idents = c("CD4+ T cells"))
DefaultAssay(cd4_t_cell) <- "RNA"

# SCTransform
features <- SplitObject(cd4_t_cell, split.by = "batch")

features <- lapply(X = features,
                   FUN = SCTransform,
                   vst.flavor = "v2",
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)

t_cell_cd4 <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)

VariableFeatures(t_cell_cd4) <- var.features
gc()

# Batch correction with Harmony
t_cell <- RunPCA(t_cell_cd4, verbose = TRUE)
t_cell <- RunHarmony(t_cell, assay.use="SCT", group.by.vars = "batch")
t_cell <- RunUMAP(t_cell, reduction = "harmony", dims = 1:30)
t_cell <- FindNeighbors(t_cell, reduction = "harmony", dims = 1:30)
t_cell <- t_cell %>% FindClusters(resolution = 0.3)

saveRDS(t_cell, "PBMC_CellRanger_6.1.2/Upload/cd4_subclusters_1.12.23.rds")

cd4_markers <- FindAllMarkers(t_cell,
                              assay = "RNA",
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              only.pos = TRUE)

write_csv(cd4_markers, "PBMC_CellRanger_6.1.2/Upload/Markers/cd4_markers_1.12.23.csv")

# CD4 T cell subsets figures
DimPlot(t_cell, label = FALSE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_t_cells_DimPLot.png", dpi = 600)

# Double-positive cells
plot_density(t_cell,
             features = c("CD8A", "CD8B"))
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subcluster_5.png", dpi = 600)

# Tregs
plot_density(t_cell,
             features = c("FOXP3", "CTLA4", "PDCD1"))
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subcluster_4.png", dpi = 600)

# Naive
plot_density(t_cell,
             features = c("CCR7", "SELL"))
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subcluster_1.png", dpi = 600)

# Cytotoxic
plot_density(t_cell,
             features = c("IL7R"))
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subcluster_2.png", dpi = 600)

plot_density(t_cell,
             features = c("NKG7", "GNLY",
                          "GZMB", "GZMA"))
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subcluster_2_extra.png", dpi = 600)




# Memory
plot_density(t_cell,
             features = c("TNFRSF4", "TIMP1"))
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subcluster_0.png", dpi = 600)

plot_density(t_cell, features = c("RORC"))

DimPlot(t_cell)

cluster.ids <- c("RORC+ cells",
                 "Naive cells I",
                 "Effector cells",
                 "Naive cells II",
                 "Tregs",
                 "CD4+ CD8+ cells",
                 "IFN cells"
)

names(cluster.ids) <- levels(t_cell)
t_cell <- RenameIdents(t_cell, cluster.ids)

DimPlot(t_cell, label = FALSE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd4_subclusters_annotated.png", dpi = 600)

saveRDS(t_cell, "PBMC_CellRanger_6.1.2/Upload/cd4_subclusters_1.12.23.rds")


######## Subcluster CD8+ T cells
rm(list=setdiff(ls(), "pbmc"))
gc()

cd8_t_cell <- subset(pbmc, idents = c("CD8+ T cells I", "CD8+ T cells II", "CD8+ T cells III"))
DefaultAssay(cd8_t_cell) <- "RNA"

# SCTransform
features <- SplitObject(cd8_t_cell, split.by = "batch")

features <- lapply(X = features,
                   FUN = SCTransform,
                   vst.flavor = "v2",
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)

t_cell_cd8 <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)

VariableFeatures(t_cell_cd8) <- var.features

gc()

# Batch correction with Harmony
t_cell <- RunPCA(t_cell_cd8, verbose = TRUE)
t_cell <- RunHarmony(t_cell, assay.use="SCT", group.by.vars = "batch")
t_cell <- RunUMAP(t_cell, reduction = "harmony", dims = 1:30)
t_cell <- FindNeighbors(t_cell, reduction = "harmony", dims = 1:30)

cd8_t_cell <- t_cell %>% FindClusters(resolution = 0.3)

cd8_markers <- FindAllMarkers(cd8_t_cell,
                              assay = "RNA",
                              logfc.threshold = 0.25,
                              min.diff.pct = 0.1,
                              test.use = "wilcox",
                              only.pos = TRUE)

write_csv(cd8_markers, "PBMC_CellRanger_6.1.2/Upload/Markers/cd8_markers_1.12.23.csv")

cluster.ids <- c("Effector cells",
                 "GZMK+ effector cells",
                 "IL7R+ TCF- cells",
                 "CD161+ effector-memory cells",
                 "Naive cells"
)

names(cluster.ids) <- levels(cd8_t_cell)
cd8_t_cell <- RenameIdents(cd8_t_cell, cluster.ids)

DimPlot(cd8_t_cell, label = FALSE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/cd8_subclusters_annotated.png", dpi = 600)

saveRDS(cd8_t_cell, "PBMC_CellRanger_6.1.2/Upload/cd8_subclusters_1.12.23.rds")

############# Subcluster NK cells
rm(list=setdiff(ls(), "pbmc"))
gc()

nk_cells <- subset(pbmc, idents = c("NK cells I", "NK cells II"))
DefaultAssay(nk_cells) <- "RNA"

# SCTransform
features <- SplitObject(nk_cells, split.by = "batch")

features <- lapply(X = features,
                   FUN = SCTransform,
                   vst.flavor = "v2",
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)

nk <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)

VariableFeatures(nk) <- var.features

gc()

# Batch correction with Harmony
subclust <- RunPCA(nk, verbose = TRUE)
subclust <- RunHarmony(subclust, assay.use="SCT", group.by.vars = "batch")
subclust <- RunUMAP(subclust, reduction = "harmony", dims = 1:30)
subclust <- FindNeighbors(subclust, reduction = "harmony", dims = 1:30)
subclust <- subclust %>% FindClusters(resolution = 0.1)

subclust %>% FeaturePlot(features = c("NKG7", "GNLY", "CD3D", "NCAM1",
                                                                         "IL2RB", "KLRC2", "GZMA", "GZMB",
                                                                         "IL7R", "SELL", "FCGR3A", "NKG2A",
                                                                         "CX3CR1", "IL18RAP", "CCL5", "CD8A"))

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/NK_subclusters.png", dpi = 600)

nk_markers <- FindAllMarkers(subclust,
                             assay = "RNA",
                             logfc.threshold = 0.25,
                             min.diff.pct = 0.2,
                             test.use = "wilcox",
                             only.pos = TRUE)
write_csv(nk_markers, "PBMC_CellRanger_6.1.2/Upload/Markers/NK_markers_1.13.23.csv")

plot_density(subclust, features = c("MKI67", "PRF1", "KLRB1",
                                    "CD38", "HAVCR2", "NCAM1"))

cluster.ids <- c("NKT cells",
                 "CD56lo NK cells I",
                 "CD56bright NK cells",
                 "CD56lo NK cells II",
                 "Proliferating NK cells")

names(cluster.ids) <- levels(subclust)
subclust <- RenameIdents(subclust, cluster.ids)

DimPlot(subclust, label = FALSE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/NK_subclusters_annotated.png", dpi = 600)

saveRDS(subclust, "PBMC_CellRanger_6.1.2/Upload/NK_subclusters_1.12.23.rds")

####### Subcluster myeloid cells
rm(list=setdiff(ls(), "pbmc"))
gc()

myeloid <- subset(pbmc, idents = c("CD14+ monocytes", "CD16+ monocytes",
                                   "CD141+ dendritic cells", "Plasmacytoid dendritic cells I",
                                   "Plasmacytoid dendritic cells II"))
DefaultAssay(myeloid) <- "RNA"

# SCTransform
features <- SplitObject(myeloid, split.by = "batch")

features <- lapply(X = features,
                   FUN = SCTransform,
                   vst.flavor = "v2",
                   method = "glmGamPoi",
                   return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = features, nfeatures = 3000)

myeloid <- merge(x = features[[1]], y = features[2:length(features)], merge.data=TRUE)

VariableFeatures(myeloid) <- var.features

gc()

# Batch correction with Harmony
subclust <- RunPCA(myeloid, verbose = TRUE)
subclust <- RunHarmony(subclust, assay.use="SCT", group.by.vars = "batch")
subclust <- RunUMAP(subclust, reduction = "harmony", dims = 1:30)
subclust <- FindNeighbors(subclust, reduction = "harmony", dims = 1:30)
subclust %>% FindClusters(resolution = 0.22) %>% FeaturePlot(features = c("FCGR3A", "CD14", "CD1C",
                                                                          "ITGAM", "LILRA4"))

ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/NK_subclusters.png", dpi = 600)

nk_markers <- FindAllMarkers(subclust,
                             assay = "RNA",
                             logfc.threshold = 0.25,
                             min.diff.pct = 0.2,
                             test.use = "wilcox",
                             only.pos = TRUE)
write_csv(nk_markers, "PBMC_CellRanger_6.1.2/Upload/Markers/NK_markers_1.13.23.csv")

plot_density(subclust, features = c("MKI67", "PRF1", "KLRB1",
                                    "CD38", "HAVCR2", "NCAM1"))

cluster.ids <- c("NKT cells",
                 "CD56lo NK cells I",
                 "CD56bright NK cells",
                 "CD56lo NK cells II",
                 "Proliferating NK cells")

names(cluster.ids) <- levels(subclust)
subclust <- RenameIdents(subclust, cluster.ids)

DimPlot(subclust, label = FALSE)
ggsave("PBMC_CellRanger_6.1.2/Upload/Figures/NK_subclusters_annotated.png", dpi = 600)

saveRDS(subclust, "PBMC_CellRanger_6.1.2/Upload/NK_subclusters_1.12.23.rds")