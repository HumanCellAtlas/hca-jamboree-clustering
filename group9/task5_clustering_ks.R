# Libraries and Code
library(Matrix)
require(MASS)
require(RANN)
require(reshape)
require(igraph)
source("ks_utils.R")

# USAGE (Given Counts_data an m x n matrix of m genes and n cells)
#cluster_ids = doGraph_clustering_KS(counts=Counts_data, num.nn = 15, num.pcs = 50)

# Load Zeisel data
zeiselData = readFormat("~/jamboree/zeisel_2015/counts/expression_mRNA_17_Aug_2014.txt")
zeiselMetaData = zeiselData[[1]]
zeiselCounts = zeiselData[[2]]

zeisel.clusters = doGraph_clustering_KS(counts=zeiselCounts, num.nn = 15, num.pcs = 50)
pdf("Zeisel_Confusion_Matrix_KS_clustering_level1.pdf",w=20,h=10)
plotConfusionMatrix(table(zeisel.clusters, zeiselMetaData$level1class))
dev.off()

pdf("Zeisel_Confusion_Matrix_KS_clustering_level2.pdf",w=40,h=10)
plotConfusionMatrix(table(zeisel.clusters, zeiselMetaData$level2class))
dev.off()

#save(list=c("zeiselCounts","zeiselMetaData"), file="ZeiselData.Rdata")

# Load PBMC data
pbmcCounts = readMM("~/jamboree/10x_pbmc/counts/from_10x/filtered_gene_bc_matrices/GRCh38/matrix.mtx")
geneids = as.character(read.csv("~/jamboree/10x_pbmc/counts/from_10x/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep="\t",header=FALSE)$V1)
cellids = as.character(read.csv("~/jamboree/10x_pbmc/counts/from_10x/filtered_gene_bc_matrices/GRCh38/barcodes.tsv", sep="\t",header=FALSE)$V1)
rownames(pbmcData) = geneids
colnames(pbmcData) = cellids

pbmc.clusters = doGraph_clustering_KS(counts=pbmcCounts, num.nn = 15, num.pcs = 50)


