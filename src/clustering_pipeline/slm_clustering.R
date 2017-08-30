# Input arguments
# arg1 --- location of the expression matrices (sparse matrix format)
# arg2 --- location of output files 
# arg3 --- name of expression matrix

args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(densityClust)
library(Matrix)
file.loc <- paste0(args[2], args[3], "_slm_results.txt")

RunSeuratClustering <- function(data.mat){
  data.mat <- t(data.mat)
  cells <- paste0("cell_", 1:ncol(data.mat))
  genes <- paste0("gene_", 1:nrow(data.mat))
  rownames(data.mat) <- genes
  colnames(data.mat) <- cells
  object <- CreateSeuratObject(raw.data = data.mat, min.cells = 3)
  object <- NormalizeData(object)
  object <- FindVariableGenes(object, do.plot = F)
  genes.use <- head(rownames(object@hvg.info), 1000)
  object <- ScaleData(object = object, genes.use = genes.use, vars.to.regress = "nUMI", model.use = "linear")
  object <- RunPCA(object, pc.genes = genes.use, pcs.compute = 50, weight.by.var = F, do.print = F)
  object <- FindClusters(object, reduction.type = "pca", dims.use = 1:20, algorithm = 3, n.start = 10)
  return(object)
}

raw.data <- readMM(args[1])
object <- RunSeuratClustering(raw.data)
results <- as.data.frame(object@ident)
results <- data.frame(results, row.names = rownames(raw.data))
write.table(results, file = file.loc, quote = F)