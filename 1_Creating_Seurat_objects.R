library(Seurat)
library(tidyverse)
library(patchwork)
library(DropletUtils)
library(glmGamPoi)
library(sctransform)
library(harmony)
library(SeuratDisk)
library(celda)

### Sample 1
counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-1/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-1/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_1_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_1_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_1_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_1_Singlet_VlnPlt.png", dpi = 300)

quantile(pbmc.singlet$nCount_RNA)
IQR(pbmc.singlet$nCount_RNA)

quantile(pbmc.singlet$nFeature_RNA, 0.75) + IQR(pbmc.singlet$nFeature_RNA)
quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_001",
      HTO_classification == "C0252" ~ "CVD_002",
      HTO_classification == "C0253" ~ "CVD_003",
      HTO_classification == "C0254" ~ "CVD_004",
      HTO_classification == "C0255" ~ "CVD_005"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_1.rds")

######### Sample 2
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-2/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-2/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_2_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_2_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_2_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_2_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_006",
      HTO_classification == "C0252" ~ "CVD_007",
      HTO_classification == "C0253" ~ "CVD_008",
      HTO_classification == "C0254" ~ "CVD_009",
      HTO_classification == "C0255" ~ "CVD_010"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_2.rds")

######### Sample 3
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-3/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-3/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_3_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_3_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_3_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_3_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_033",
      HTO_classification == "C0252" ~ "CVD_034",
      HTO_classification == "C0253" ~ "CVD_035",
      HTO_classification == "C0254" ~ "CVD_036",
      HTO_classification == "C0255" ~ "CVD_040"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_3.rds")

######### Sample 4
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-4/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-4/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_4_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_4_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_4_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_4_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_038",
      HTO_classification == "C0252" ~ "CVD_039",
      HTO_classification == "C0253" ~ "CVD_011",
      HTO_classification == "C0254" ~ "CVD_012",
      HTO_classification == "C0255" ~ "CVD_013"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_4.rds")

######### Sample 5
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-5/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-5/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_5_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_5_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_5_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_5_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_041",
      HTO_classification == "C0252" ~ "CVD_023",
      HTO_classification == "C0253" ~ "CVD_024",
      HTO_classification == "C0254" ~ "CVD_025",
      HTO_classification == "C0255" ~ "CVD_021"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_5.rds")

######### Sample 6
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-6/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-6/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_6_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_6_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_6_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_6_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_031",
      HTO_classification == "C0252" ~ "CVD_032",
      HTO_classification == "C0253" ~ "CVD_017",
      HTO_classification == "C0254" ~ "CVD_018",
      HTO_classification == "C0255" ~ "CVD_019"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_6.rds")

######### Sample 7
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-7/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-7/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_7_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_7_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_7_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_7_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_026",
      HTO_classification == "C0252" ~ "CVD_027",
      HTO_classification == "C0253" ~ "CVD_028",
      HTO_classification == "C0254" ~ "CVD_029",
      HTO_classification == "C0255" ~ "CVD_030"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_7.rds")

######### Sample 8
rm(list = ls())
gc()

counts <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-8/filtered_feature_bc_matrix.h5")
counts.raw <- Read10X_h5("PBMC_CellRanger_6.1.2/8256-SB-8/raw_feature_bc_matrix.h5")

# Correcting for ambient RNA using DecontX
sce <- SingleCellExperiment(list(counts = counts$`Gene Expression`))
sce.raw <- SingleCellExperiment(list(counts = counts.raw$`Gene Expression`))

sce <- decontX(sce, background = sce.raw)

# Merging with HTOs
pbmc.hashtag <- CreateSeuratObject(counts = round(decontXcounts(sce)), min.cells = 3, min.features = 100)
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_8_Full_VlnPlt.png", dpi = 300)

plot1 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_8_nFeatures_nCount.png", dpi = 300, height = 7, width = 17)

pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = counts.raw$`Antibody Capture`[,colnames(pbmc.hashtag)])
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]]), ncol = 3)
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_8_HTO_ridgeplot.png", dpi = 300)

Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

VlnPlot(pbmc.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave("PBMC_CellRanger_6.1.2/Upload/QC/Sample_8_Singlet_VlnPlt.png", dpi = 300)

x <- (quantile(pbmc.singlet$nFeature_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nFeature_RNA))
y <- (quantile(pbmc.singlet$nCount_RNA, 0.75) + 1.5*IQR(pbmc.singlet$nCount_RNA))

pbmc.singlet <- subset(pbmc.singlet,
                       subset = nFeature_RNA > 100 & nFeature_RNA < x & nCount_RNA < y & percent.mt < 10)

all(pbmc.singlet@meta.data$HTO_classification == pbmc.singlet@meta.data$hash.ID)

pbmc.singlet@meta.data <- pbmc.singlet@meta.data %>%
  dplyr::mutate(
    Sample_ID = case_when(
      HTO_classification == "C0251" ~ "CVD_014",
      HTO_classification == "C0252" ~ "CVD_015",
      HTO_classification == "C0253" ~ "CVD_016",
      HTO_classification == "C0254" ~ "CVD_020",
      HTO_classification == "C0255" ~ "CVD_022"
    ))

saveRDS(pbmc.singlet, "PBMC_CellRanger_6.1.2/Upload/decontX_sample_8.rds")