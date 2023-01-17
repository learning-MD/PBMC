library(edgeR)
library(ggplot2)
library(ggpubr)
library(testit)
library(tools)
library(Matrix.utils)
library(limma)

# Load data
counts <- read_csv("PBMC_CellRanger_6.1.2/Upload/Pseudocounts/Platelets.pseudo_count.csv")
counts <- counts %>% dplyr::rename(gene = `...1`)

# Remove ribosomal and mitochondrial genes
counts <- counts[!grepl("^RP[SL]", counts$gene),]
counts <- counts[!grepl("MT-", counts$gene),]

metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
metadata$prednisone <- as.factor(metadata$prednisone)

# Limma
# Comparing CAV vs non-CAV
y <- DGEList(counts = counts[,2:41], genes = counts$gene, group = metadata$prednisone)
y$samples
y

keep.exprs <- filterByExpr(y, group=metadata$prednisone)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

y <- calcNormFactors(y)
y$samples

# Evaluating for gender differences
Gender <- substring(metadata$sex,1,1)
plotMDS(y, labels=Gender, top=50, col=ifelse(Gender=="m","blue","red"), gene.selection="common", prior.count = 5)

# Designing matrix
mm <- model.matrix(~0 + prednisone+age+sex, data = metadata)
mm

# Voom
v <- voom(y, mm, plot = TRUE, normalize = "quantile")

fit <- lmFit(v, mm)  

contr <- makeContrasts(prednisone1 - prednisone0, levels = colnames(mm))

tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)

tt <- topTable(tmp, n = Inf, adjust.method = "BH")

df <- data.frame(gene = tt$genes,
                 logFC = tt$logFC,
                 pval = tt$P.Value,
                 padf = tt$adj.P.Val
)


write_csv(df, "PBMC_CellRanger_6.1.2/Upload/DGE/Prednisone_v_nonPred/Platelets.csv")
rm(list=setdiff(ls(), c("pbmc")))
gc()

#############################################################################################
# Modifications made for when a sample does NOT contain a cluster
# Remove that sample from metadata
metadata <- read_csv("PBMC_CellRanger_6.1.2/Upload/pbmc_metadata.csv")
metadata <- metadata %>% filter(Sample_ID != "CVD_010" & Sample_ID != "CVD_019")

# Limma
# Comparing CAV vs non-CAV
y <- DGEList(counts = counts[,2:39], genes = counts$gene, group = metadata$prednisone)
y$samples
y
