# Create pseudocounts for each cluster
sumcount<-function(ct_count, groupings){
  rescount<-t(aggregate.Matrix(t(ct_count), groupings=groupings, fun="sum"))
  psum<-rowSums(rescount)
  rescount<-rescount[psum>0,]
  return(rescount)
}

pbmc <- readRDS("PBMC_CellRanger_6.1.2/Upload/pbmc_annotated_exclude_TCR.rds")
pbmc@meta.data$cluster_name <- Idents(pbmc)

clusterDf<-pbmc@meta.data

cluster_name = "cluster_name"
cts = as.character(unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T),
                                    cluster_name]))
prefixList<-gsub('[ /:_]+', '_', cts)

res_files=c()

idx<-1
for (idx in c(1:length(cts))){
  ct = cts[idx]
  cat(ct, "\n")
  
  prefix = prefixList[idx]
  
  clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
  de_obj<-subset(pbmc, cells=rownames(clusterCt))
  
  ct_count<-de_obj@assays$RNA@counts
  groupings<-unlist(de_obj$Sample_ID)
  p_count<-sumcount(ct_count, groupings)
  
  p_file=paste0(prefix, ".pseudo_count.csv")
  
  write.csv(p_count, p_file)
  
  res_files<-c(res_files, file_path_as_absolute(p_file))
}
