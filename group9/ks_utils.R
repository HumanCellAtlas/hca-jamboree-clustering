# Function for estimating variable genes and graph clustering
require(MASS)
require(RANN)
require(reshape)
require(igraph)

new.meanCVfit <- function(count.data, diffCV.cutoff=3, mean.min=.005, mean.max=1000, main.use="", do.plot=T, num.sd=1.5) {
  
  # Format for the plot
  if (do.plot) par(mfrow=c(1,3), oma=c(0,0,2,0))
  
  # Calculate empirical mean, var and CV
  mean_emp <- apply(count.data, 1, mean)
  var_emp <- apply(count.data, 1, var)
  cv_emp <- sqrt(var_emp) / mean_emp
  
  # Fit gamma distribution to 'size factors' to build negative binomial
  a <- colSums(count.data)
  size_factor <- a/mean(a)
  fit <- fitdistr(size_factor, "Gamma")
  if (do.plot) {
    print(fit)
    hist(size_factor, 50, probability=TRUE, xlab="UMIs per cell / mean UMIs per cell", main=paste0("Size Factors & Gamma Fit (a=", round(fit$estimate[1], 2), ")"))
    curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.999), add=TRUE, col="red")
  }
  
  # Create null negative binomial distribution
  #   Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
  #   then cX ~ Gamma(a, b/c)
  a_i <- rep(fit$estimate[1], length(mean_emp))
  names(a_i) <- names(mean_emp)
  b_i <- fit$estimate[2] / mean_emp
  names(b_i) <- names(mean_emp)
  mean_NB <- a_i / b_i
  var_NB <- a_i*(1+b_i) / (b_i^2)
  cv_NB <- sqrt(var_NB)/mean_NB
  
  # Calculate difference in genes' CV and null CV
  diffCV = log(cv_emp) - log(cv_NB)
  if (!is.null(num.sd)){
    diffCV.cutoff = mean(diffCV) + num.sd*sd(diffCV)
  }
  

  if (do.plot) {
    hist(diffCV, 150, xlab="log(CVgene / CVnull)", main="Diff CV")
    abline(v=diffCV.cutoff, lty=2, col='red')
  }
  
 
  pass.cutoff <- names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > mean.min & mean_emp < mean.max))]
  
  # Plot variable gene selection
  if (do.plot) {
    plot(mean_emp,cv_emp,pch=16,cex=0.5,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy", main = "Selection of Variable Genes")
    curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
    or = order(mean_NB)
    lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
    points(mean_emp[pass.cutoff], cv_emp[pass.cutoff], pch=16, cex=0.5, col='blue')
    title(main.use, outer=T)
  }
  
  # Return list of variable genes
  return(pass.cutoff)
}

readFormat <- function(infile) { 
  metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] # First column is empty.
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]
  metadata <- as.data.frame(t(metadata))
  counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] # First column after row names is some useless filler.
  counts <- as.matrix(counts)
  return(list(metadata=metadata, counts=counts))
}


# Graph clustering

doGraph_clustering_KS =  function(counts=NULL, num.nn = 30, min.cells = 30, min.genes = 400, num.sd = 1.5, num.pcs=25){

  
  counts_filt = counts[rowSums(counts > 0) > min.cells,] 
  print(paste0("Removing ", nrow(counts) - nrow(counts_filt)," genes that are expressed in fewer than ", min.cells, " cells"))
  counts_filt = counts_filt[,colSums(counts > 0) > min.genes]
  print(paste0("Removing ", ncol(counts) - ncol(counts_filt)," cells that express fewer than ", min.genes, " genes"))
  
  # perform median normalization
  print("Performing median-log normalization. Convering raw counts to log(TPX+1)")
  med_trans = median(colSums(counts_filt))
  normData = scale(counts_filt, center=FALSE, scale=colSums(counts_filt))
  logNormCounts = log(med_trans*normData + 1)
  
  var_genes = new.meanCVfit(counts_filt, mean.min = 0.005, num.sd = num.sd)
  
  # Perform PCA
  print(paste0("Performing PCA ..."))
  pca.out = princomp(t(logNormCounts[var_genes,]))
  
  print(paste0("Using ", num.pcs, " PCs for clustering. Change num.pcs parameter to the desired value" ))
  
  print(paste0("Computing an kNN graph with k=", num.nn,", and jaccard weighting."))
  Adj = getAdjMatrix(pca.out$scores[,c(1:num.pcs)],nn=num.nn,edge.weights=FALSE,do.jaccard=TRUE,do.sparse=TRUE)
  print("Performing Louvain clustering ")
  g <- graph.adjacency(Adj, weighted=TRUE, mode="undirected")
  graph.out = cluster_louvain(g)
  clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
  names(clust.assign) = graph.out$names
  # Reassign cluster labels to order clusters in decreasing numbers of cells
  k=order(table(clust.assign), decreasing = TRUE)
  new.levels = rep(1,length(unique(graph.out$membership)))
  new.levels[k] = 1:length(unique(graph.out$membership))
  levels(clust.assign) = new.levels
  clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
  return(clust.assign)   
}

  
getAdjMatrix=function(X,nn=30,edge.weights=FALSE,do.jaccard=TRUE,do.sparse=TRUE) {
    if (dim(X)[1] < 500){
      searchtype.use="standard"
      treetype.use = "kd"
    } else {
      searchtype.use="priority"
      treetype.use = "bd"
    }
    nearest=nn2(X,X,k=nn+1, treetype = treetype.use, searchtype=searchtype.use)
    print("Found nearest neighbors")
    nearest$nn.idx = nearest$nn.idx[,-1]
    nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
    
    if (edge.weights){
      nearest$nn.dists = nearest$nn.dists / median(nearest$nn.dists)
      nearest$nn.sim = 1 / nearest$nn.dists  #Convert to a similarity score
    } else { 
      nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
    }
    
    edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
    edges$B = edges$C; edges$C=1
    
    #Remove repetitions
    edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
    

    if (do.jaccard){
        NN = nearest$nn.idx
        jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
        
        edges$C = jaccard_dist
        edges = subset(edges, C != 0)
        edges$C = edges$C/max(edges$C)
      } 
      
      Adj = Matrix(0, nrow=nrow(X), ncol=nrow(X), sparse=do.sparse)
      Adj[cbind(edges$A,edges$B)] = edges$C
      Adj[cbind(edges$B,edges$A)] = edges$C
      rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
      return(Adj)
    
    
}

# Plot confusion matrix
plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, cols.use=gray.colors(10), max.size=5, ylab.use="Known", xlab.use="Predicted"){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  #X = X[rev(1:dim(X)[1]),]
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  X$Predicted = as.factor(X$Predicted)
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low ="blue",   high = "red", limits=c(0, 100 ))+scale_size(range = c(1, max.size))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))  
  print(p)
}