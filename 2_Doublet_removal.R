library(Seurat)
library(scDblFinder)
library(Matrix)
library(tidyverse)

# Read in metadata
meta <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")

# Sample 1
sample_1 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_1.rds")
sample_1_sce <- as.SingleCellExperiment(sample_1)

sce_1 <- scDblFinder(sample_1_sce, samples = "Sample_ID")
table(sce_1$scDblFinder.class, sce_1$scDblFinder.sample)

sce_1_seurat <- as.Seurat(sce_1)
DefaultAssay(sce_1_seurat) <- "RNA"

sce_1_seurat <- subset(sce_1_seurat, subset = scDblFinder.class == "singlet")
sce_1_seurat@meta.data <- left_join(sce_1_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 2
sample_2 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_2.rds")
sample_2_sce <- as.SingleCellExperiment(sample_2)

sce_2 <- scDblFinder(sample_2_sce, samples = "Sample_ID")
table(sce_2$scDblFinder.class, sce_2$scDblFinder.sample)

sce_2_seurat <- as.Seurat(sce_2)
DefaultAssay(sce_2_seurat) <- "RNA"

sce_2_seurat <- subset(sce_2_seurat, subset = scDblFinder.class == "singlet")
sce_2_seurat@meta.data <- left_join(sce_2_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 3
sample_3 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_3.rds")
sample_3_sce <- as.SingleCellExperiment(sample_3)

sce_3 <- scDblFinder(sample_3_sce, samples = "Sample_ID")
table(sce_3$scDblFinder.class, sce_3$scDblFinder.sample)

sce_3_seurat <- as.Seurat(sce_3)
DefaultAssay(sce_3_seurat) <- "RNA"

sce_3_seurat <- subset(sce_3_seurat, subset = scDblFinder.class == "singlet")
sce_3_seurat@meta.data <- left_join(sce_3_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 4
sample_4 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_4.rds")
sample_4_sce <- as.SingleCellExperiment(sample_4)

sce_4 <- scDblFinder(sample_4_sce, samples = "Sample_ID")
table(sce_4$scDblFinder.class, sce_4$scDblFinder.sample)

sce_4_seurat <- as.Seurat(sce_4)
DefaultAssay(sce_4_seurat) <- "RNA"

sce_4_seurat <- subset(sce_4_seurat, subset = scDblFinder.class == "singlet")
sce_4_seurat@meta.data <- left_join(sce_4_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 5
sample_5 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_5.rds")
sample_5_sce <- as.SingleCellExperiment(sample_5)

sce_5 <- scDblFinder(sample_5_sce, samples = "Sample_ID")
table(sce_5$scDblFinder.class, sce_5$scDblFinder.sample)

sce_5_seurat <- as.Seurat(sce_5)
DefaultAssay(sce_5_seurat) <- "RNA"

sce_5_seurat <- subset(sce_5_seurat, subset = scDblFinder.class == "singlet")
sce_5_seurat@meta.data <- left_join(sce_5_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 6
sample_6 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_6.rds")
sample_6_sce <- as.SingleCellExperiment(sample_6)

sce_6 <- scDblFinder(sample_6_sce, samples = "Sample_ID")
table(sce_6$scDblFinder.class, sce_6$scDblFinder.sample)

sce_6_seurat <- as.Seurat(sce_6)
DefaultAssay(sce_6_seurat) <- "RNA"

sce_6_seurat <- subset(sce_6_seurat, subset = scDblFinder.class == "singlet")
sce_6_seurat@meta.data <- left_join(sce_6_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 7
sample_7 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_7.rds")
sample_7_sce <- as.SingleCellExperiment(sample_7)

sce_7 <- scDblFinder(sample_7_sce, samples = "Sample_ID")
table(sce_7$scDblFinder.class, sce_7$scDblFinder.sample)

sce_7_seurat <- as.Seurat(sce_7)
DefaultAssay(sce_7_seurat) <- "RNA"

sce_7_seurat <- subset(sce_7_seurat, subset = scDblFinder.class == "singlet")
sce_7_seurat@meta.data <- left_join(sce_7_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Sample 8
sample_8 <- read_rds("PBMC_CellRanger_6.1.2/Upload/decontX_sample_8.rds")
sample_8_sce <- as.SingleCellExperiment(sample_8)

sce_8 <- scDblFinder(sample_8_sce, samples = "Sample_ID")
table(sce_8$scDblFinder.class, sce_8$scDblFinder.sample)

sce_8_seurat <- as.Seurat(sce_8)
DefaultAssay(sce_8_seurat) <- "RNA"

sce_8_seurat <- subset(sce_8_seurat, subset = scDblFinder.class == "singlet")
sce_8_seurat@meta.data <- left_join(sce_8_seurat@meta.data,
                                    meta,
                                    by = "Sample_ID")

# Adding batch
sce_1_seurat@meta.data$batch <- "batch_1"

sce_2_seurat@meta.data$batch <- "batch_2"

sce_3_seurat@meta.data$batch <- "batch_3"

sce_4_seurat@meta.data$batch <- "batch_4"

sce_5_seurat@meta.data$batch <- "batch_5"

sce_6_seurat@meta.data$batch <- "batch_6"

sce_7_seurat@meta.data$batch <- "batch_7"

sce_8_seurat@meta.data$batch <- "batch_8"

# All 8 merged samples
pbmc.combined <- merge(sce_1_seurat, y =
                         c(sce_2_seurat,
                           sce_3_seurat,
                           sce_4_seurat,
                           sce_5_seurat,
                           sce_6_seurat,
                           sce_7_seurat,
                           sce_8_seurat))

saveRDS(pbmc.combined, "PBMC_CellRanger_6.1.2/Upload/decontX_all_combined.rds")
