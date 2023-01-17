library(edgeR)
library(ggplot2)
library(ggpubr)
library(testit)
library(tools)
library(Matrix.utils)
library(limma)
library(tidyverse)

#### CD8 T cells
# Set folder path
folder_path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD8/"

# Get list of csv files in folder
csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Loop through csv files
for (csv_file in csv_files) {
  
  # Get file name without path and extension
  file_name <- gsub("^.*/|\\.pseudo_count.csv$", "", csv_file)
  
  # Load data
  counts <- read_csv(csv_file)
  counts <- counts %>% dplyr::rename(gene = `...1`)
  
  # Remove ribosomal and mitochondrial genes
  counts <- counts[!grepl("^RP[SL]", counts$gene),]
  counts <- counts[!grepl("MT-", counts$gene),]
  
  # Load metadata
  metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
  
  # Pare metadata samples to the ones represented inthe counts column
  metadata <- metadata %>% filter(Sample_ID %in% colnames(counts))
  
  # Limma
  # Comparing CAV vs non-CAV
  y <- DGEList(counts = dplyr::select(counts, -gene), genes = counts$gene, group = metadata$group)
  
  keep.exprs <- filterByExpr(y, group=metadata$group)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  # Designing matrix
  mm <- model.matrix(~0 + group+age+sex, data = metadata)
  
  # Voom
  v <- voom(y, mm, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, mm)  
  contr <- makeContrasts(groupcav - groupcontrol, levels = colnames(mm))
  tmp <- contrasts.fit(fit, contrasts = contr)
  tmp <- eBayes(tmp)
  tt <- topTable(tmp, n = Inf, adjust.method = "BH")
  df <- data.frame(gene = tt$genes,
                   logFC = tt$logFC,
                   pval = tt$P.Value,
                   padf = tt$adj.P.Val
  )
  
  # Save new csv file
  write_csv(df, paste0("PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD8/CAV_group/", file_name, ".csv"))
}


####################### NK cells
# Set folder path
folder_path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/NK_cells/"

# Get list of csv files in folder
csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Loop through csv files
for (csv_file in csv_files) {
  
  # Get file name without path and extension
  file_name <- gsub("^.*/|\\.pseudo_count.csv$", "", csv_file)
  
  # Load data
  counts <- read_csv(csv_file)
  counts <- counts %>% dplyr::rename(gene = `...1`)
  
  # Remove ribosomal and mitochondrial genes
  counts <- counts[!grepl("^RP[SL]", counts$gene),]
  counts <- counts[!grepl("MT-", counts$gene),]
  
  # Load metadata
  metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
  metadata$group <- as.factor(metadata$group)
  
  # Pare metadata samples to the ones represented inthe counts column
  metadata <- metadata %>% filter(Sample_ID %in% colnames(counts))
  
  # Limma
  # Comparing CAV vs non-CAV
  y <- DGEList(counts = dplyr::select(counts, -gene), genes = counts$gene, group = metadata$group)
  
  keep.exprs <- filterByExpr(y, group=metadata$group)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  # Designing matrix
  mm <- model.matrix(~0 + group+age+sex, data = metadata)
  
  # Voom
  v <- voom(y, mm, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, mm)  
  contr <- makeContrasts(groupcav - groupcontrol, levels = colnames(mm))
  tmp <- contrasts.fit(fit, contrasts = contr)
  tmp <- eBayes(tmp)
  tt <- topTable(tmp, n = Inf, adjust.method = "BH")
  df <- data.frame(gene = tt$genes,
                   logFC = tt$logFC,
                   pval = tt$P.Value,
                   padf = tt$adj.P.Val
  )
  
  # Save new csv file
  write_csv(df, paste0("PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/NK_cells/CAV_group/", file_name, ".csv"))
}

################## CD4 T cells
# Set folder path
folder_path <- "PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD4/"

# Get list of csv files in folder
csv_files <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Loop through csv files
for (csv_file in csv_files) {
  
  # Get file name without path and extension
  file_name <- gsub("^.*/|\\.pseudo_count.csv$", "", csv_file)
  
  # Load data
  counts <- read_csv(csv_file)
  counts <- counts %>% dplyr::rename(gene = `...1`)
  
  # Remove ribosomal and mitochondrial genes
  counts <- counts[!grepl("^RP[SL]", counts$gene),]
  counts <- counts[!grepl("MT-", counts$gene),]
  
  # Load metadata
  metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
  metadata$group <- as.factor(metadata$group)

  # Pare metadata samples to the ones represented inthe counts column
  metadata <- metadata %>% filter(Sample_ID %in% colnames(counts))
  
  # Limma
  # Comparing CAV vs non-CAV
  y <- DGEList(counts = dplyr::select(counts, -gene), genes = counts$gene, group = metadata$group)
  
  keep.exprs <- filterByExpr(y, group=metadata$group)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  # Designing matrix
  mm <- model.matrix(~0 + group+age+sex, data = metadata)
  
  # Voom
  v <- voom(y, mm, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, mm)  
  contr <- makeContrasts(groupcav - groupcontrol, levels = colnames(mm))
  tmp <- contrasts.fit(fit, contrasts = contr)
  tmp <- eBayes(tmp)
  tt <- topTable(tmp, n = Inf, adjust.method = "BH")
  df <- data.frame(gene = tt$genes,
                   logFC = tt$logFC,
                   pval = tt$P.Value,
                   padf = tt$adj.P.Val
  )
  
  # Save new csv file
  write_csv(df, paste0("PBMC_CellRanger_6.1.2/Upload/Subclusters_DGE/CD4/CAV_group/", file_name, ".csv"))
}
