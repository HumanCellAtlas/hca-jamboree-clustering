## Questions? -> nf@ebi.ac.uk

install.packages <- function() {
    source("https://bioconductor.org/biocLite.R")
    biocLite("SC3",suppressUpdates=TRUE)
    biocLite("scater",suppressUpdates=TRUE)
    biocLite("agricolae",suppressUpdates=TRUE)
}

DEBUG<-FALSE

## imatrix - count matrix
## ks - number of clusters
gen.clusters.sc3 <- function(imatrix,ks,biology=FALSE,...) {
    
    library(scater)
    library(SC3)

    stopifnot(is.numeric(ks))
    tmp <- imatrix
    ## SCESEt object
    ##sceset <- newSCESet(countData = tmp, phenoData = NULL, logExprsOffset = 1)
    sceset <- newSCESet(countData = tmp)
    ## QC 
    sceset <- calculateQCMetrics(sceset)

    if (DEBUG)
        plotPCA(sceset)

    ## 
    sceset <- sc3(sceset, ks = ks, biology = biology, ...)

    sc3_plot_consensus(sceset,k=max(ks))
    if(DEBUG)
        sc3_plot_consensus(sceset, k = max(ks))
    sc3_plot_cluster_stability(sceset,k=max(ks))
    p_data <- pData(sceset)
    label <- paste0("sc3_",max(ks),"_clusters")
    clus <- p_data[,label]
    names(clus) <- rownames(p_data)
    return(clus)
}

gen.clusters.hclust <- function(imatrix,cor.method="spearman",bootstrap=100,distance="euclidean",...) {
    library("agricolae")

    ## assumes raw counts matrix
    imatrix.cor <- round(cor(imatrix,method=cor.method),2)
    ##d <- dist(imatrix.cor,"euclidean")
    con <- consensus(imatrix.cor,distance=distance,nboot=bootstrap,...)
    ##
    ## TODO: pick a better h
    h <- max(con$dendrogram$height)-max(con$dendrogram$height)/3
    x <- cutree(con$dendrogram,h=h)
    return(x)
}
## Run
## c1 <- gen.clusters.sc3(imatrix,ks=9,n_cores=3)
## c2 <- gen.clusters.hclust(imatrix,cor.method="spearman",nboot=100)

