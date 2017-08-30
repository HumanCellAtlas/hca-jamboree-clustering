args <- commandArgs(trailingOnly = TRUE)
library(monocle)
data.mat <- args[1] 
file.loc <- paste0(args[2], args[3], "_tsne_results.txt")

tsne_monocle <- function(mat) {
  mat <- t(mat)
  row.names(mat) <- paste0("gene_", 1:nrow(mat))
  colnames(mat) <- paste0("cell_", 1:ncol(mat))
  dhsinfo <- data.frame(x = 1:nrow(mat))
  row.names(dhsinfo) <- row.names(mat)
  
  cellinfo <- data.frame(x = 1:ncol(mat))
  row.names(cellinfo) <- colnames(mat)
  
  fd <- new("AnnotatedDataFrame", data = dhsinfo)
  pd <- new("AnnotatedDataFrame", data = cellinfo)
  GM_cds <-  newCellDataSet(as(mat, "sparseMatrix"),
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit=1)
  GM_cds <- detectGenes(GM_cds, min_expr = 0.1)
  GM_cds <- estimateSizeFactors(GM_cds)
  GM_cds <- estimateDispersions(GM_cds)
  
  disp_table <- dispersionTable(GM_cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  GM_cds <- setOrderingFilter(GM_cds, unsup_clustering_genes$gene_id)
  
  
  GM_cds <- reduceDimension(GM_cds, max_components = 2, num_dim = 5,
                            reduction_method = 'tSNE', verbose = T)
  GM_cds <- clusterCells(GM_cds, method="densityPeak")
  return(pData(GM_cds)$Cluster)
}

data.mat <- readMM(args[1])
results <- tsne_monocle(data.mat)
results <- data.frame(results, row.names = rownames(data.mat))
write.table(results, file = file.loc, quote = F)