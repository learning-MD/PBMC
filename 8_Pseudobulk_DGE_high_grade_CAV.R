library(edgeR)
library(ggplot2)
library(ggpubr)
library(testit)
library(tools)
library(Matrix.utils)
library(limma)

# Adjusting for prednisone and removing all CAV-1

# Set folder path
folder_path <- "PBMC_CellRanger_6.1.2/Upload/Pseudocounts"

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
  
  # Load metadata and include only CAV-0, CAV-2, and CAV-3
  metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
  metadata <- metadata %>% filter(
    cav_grade == 0 | cav_grade == 2 | cav_grade == 3
  )
  metadata$high_grade_cav <- as.factor(metadata$high_grade_cav)
  metadata$group <- as.factor(metadata$group)
  metadata$prednisone <- as.factor(metadata$prednisone)
  
  # Remove columns in the count matrix associated with CAV-1
  counts <- counts %>% select(gene, any_of(metadata$Sample_ID))
  
  # Pare metadata samples to the ones represented inthe counts column
  metadata <- metadata %>% filter(Sample_ID %in% colnames(counts))
  
  # Limma
  # Comparing CAV vs non-CAV
  y <- DGEList(counts = dplyr::select(counts, -gene), genes = counts$gene, group = metadata$group)
  
  keep.exprs <- filterByExpr(y, group=metadata$group)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  # Designing matrix
  mm <- model.matrix(~0 + group+age+sex+prednisone, data = metadata)
  
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
  write_csv(df, paste0("PBMC_CellRanger_6.1.2/Upload/DGE/High_grade_vs_nonCAV/Pred_adjusted/", file_name, ".csv"))
}
