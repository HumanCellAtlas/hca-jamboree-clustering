# clustering functions: 
# input should be cm - an integer matrix of counts; mc.cores - an integer defining how many cores to use
# 1) group9_cluster_CH1(cm, mc.cores=1) 
# 2) group9_cluster_CH2(cm, mc.cores=1)
# 3) group9_cluster_CH3(cm)
# 4) group9_cluster_CH4(cm) 

library('Matrix')
library('MASS')
library('diffusionMap')
library('FNN')
library('igraph')
library('onlinePCA')
library('parallel')

theta.reg <- function(cm, regressors, min.theta=0.01, bins=12, reg.fac=1, bandwidth=0.05) {
  b.id <- (1:nrow(cm)) %% max(1, bins, na.rm=TRUE) + 1
  cat(sprintf('get regularized theta estimate for %d genes and %d cells\n', nrow(cm), ncol(cm)))
  cat(sprintf('processing %d bins with ca %d genes in each\n', bins, round(nrow(cm)/bins, 0)))
  theta.estimate <- rep(NA, nrow(cm))
  for (bin in sort(unique(b.id))) {
    sel.g <- which(b.id == bin)
    bin.theta.estimate <- unlist(mclapply(sel.g, function(i) {
      as.numeric(theta.ml(cm[i, ], glm(cm[i, ] ~ ., data = regressors, family=poisson)$fitted))
    }), use.names = FALSE)
    theta.estimate[sel.g] <- bin.theta.estimate
    cat(sprintf('%d ', bin))
  }
  cat('done\n')
  raw.mean <- apply(cm, 1, mean)
  log.raw.mean <- log10(raw.mean)
  var.estimate <- raw.mean + raw.mean^2/theta.estimate
  
  fit <- loess(log10(var.estimate) ~ log.raw.mean, span=0.33)
  #points(log.raw.mean, fit$fitted, col='red', pch='.')
  theta.fit <- raw.mean^2 / (10^fit$fitted - raw.mean)
  #var.fit <- raw.mean + raw.mean^2/theta.fit
  
  to.fix <- theta.fit <= min.theta | is.infinite(theta.fit)
  if (any(to.fix)) {
    cat('Fitted theta below', min.theta, 'for', sum(to.fix), 'genes, setting them to', min.theta, '\n')
    theta.fit[to.fix] <- min.theta
  }
  names(theta.fit) <- rownames(cm)
  return(theta.fit)
}

nb.residuals.glm <- function(y, regression.mat, fitted.theta, gene) {
  #return(1:length(y))
  fit <- 0
  try(fit <- glm(y ~ ., data = regression.mat, family=negative.binomial(theta=fitted.theta)), silent=TRUE)
  if (class(fit)[1] == 'numeric') {
    message(sprintf('glm and family=negative.binomial(theta=%f) failed for gene %s; falling back to scale(log10(y+1))', 
                    fitted.theta, gene))
    return(scale(log10(y+1))[, 1])
  }
  return(residuals(fit, type='pearson'))
}

norm.nb.reg <- function(cm, regressors, min.theta=0.01, bins=64, theta.fit=NA, pr.th=NA, save.theta.fit=c()) {
  cat('Normalizing data using regularized NB regression\n')
  cat('explanatory variables:', colnames(regressors), '\n')
  if (any(is.na(theta.fit))) {
    theta.fit <- theta.reg(cm, regressors, min.theta, bins)
    if (is.character(save.theta.fit)) {
      save(theta.fit, file=save.theta.fit)
    }
  }
  b.id <- (1:nrow(cm)) %% max(1, bins, na.rm=TRUE) + 1
  cat('Running NB regression\n')
  res <- matrix(NA, nrow(cm), ncol(cm), dimnames=dimnames(cm))
  for (bin in sort(unique(b.id))) {
    sel.g <- rownames(cm)[b.id == bin]
    expr.lst <- mclapply(sel.g, function(gene) nb.residuals.glm(cm[gene, ], regressors, theta.fit[gene], gene), mc.preschedule = FALSE)
    res[sel.g, ] <- do.call(rbind, expr.lst)
    cat(sprintf('%d ', bin))
  }
  cat('done\n')
  if (!any(is.na(pr.th))) {
    res[res > pr.th] <- pr.th
    res[res < -pr.th] <- -pr.th
  }
  attr(res, 'theta.fit') <- theta.fit
  return(res)
}

jacc.sim <- function(mat, k) {
  # generate a sparse nearest neighbor matrix
  #require('FNN')
  nn.indices <- get.knn(mat, k)$nn.index
  j <- as.numeric(t(nn.indices))
  i <- ((1:length(j))-1) %/% k + 1
  nn.mat <- sparseMatrix(i=i, j=j, x=1)
  rm(nn.indices, i, j)
  # turn nn matrix into SNN matrix and then into Jaccard similarity
  snn <- nn.mat %*% t(nn.mat)
  snn.summary <- summary(snn)
  snn <- sparseMatrix(i=snn.summary$i, j=snn.summary$j, x=snn.summary$x/(2*k-snn.summary$x))
  rm(snn.summary)
  return(snn)
}


cluster.the.data.simple <- function(cm, expr, k, sel.g=NA, d.meth='cor', do.scale=TRUE, do.find.markers=FALSE, 
                                    min.mean=0.001, min.cells=3, z.th=1, ev.red.th=0.04, in.evdim=NA, n.perms=1, seed=NULL, 
                                    max.dim=50) {
  goi <- rownames(expr)[apply(cm[rownames(expr), ]>0, 1, sum) >= min.cells & apply(cm[rownames(expr), ], 1, mean) >= min.mean]
  sspr <- apply(expr[goi, ]^2, 1, sum)
  sel.g <- goi[scale(sqrt(sspr)) > z.th]

  cat(sprintf('Selected %d variable genes\n', length(sel.g)))
  sel.g <- intersect(sel.g, rownames(expr))
  cat(sprintf('%d of these are in expression matrix.\n', length(sel.g)))
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  if (sum(is.na(expr[sel.g, ])) > 0) {
    dmat <- 1 - cor(expr[sel.g, ], use = 'pairwise.complete.obs')
  } else {
    dmat <- 1 - cor(expr[sel.g, ])
  }

  max.dim <- min(max.dim, nrow(dmat)/2)
  dmap <- diffuse(dmat, neigen=max.dim, maxdim=max.dim)
  dm <- dmap$X
  ev <- dmap$eigenvals
  
  #ev.red <- diff(ev)/-ev[-length(ev)]
  ev.red <- ev/sum(ev)
  #evdim <- rev(which(scale(ev.red) > ev.red.th))[1]
  evdim <- rev(which(ev.red > ev.red.th))[1]
  
  plot(ev, ylim=c(0, max(ev)))
  abline(v=evdim + 0.5, col='blue')
  
  cat('Using', evdim, 'significant DM coordinates\n')
  evdim <- max(2, evdim, na.rm=TRUE)
  
  # for now just use the first random matrix as background
  if (do.scale) {
    dm <- scale(dm[, 1:evdim])
  } else {
    dm <- dm[, 1:evdim]
  }
  sim.mat <- jacc.sim(dm, k)
  
  gr <- graph_from_adjacency_matrix(sim.mat, mode='undirected', weighted=TRUE, diag=FALSE)
  gr$layout <- layout.fruchterman.reingold(gr)
  cl <- as.numeric(membership(cluster_louvain(gr)))
  
  results <- list()
  results$dm <- dm
  results$dm.orig <- dmap
  results$clustering <- cl
  results$sel.g <- sel.g
  results$sim.mat <- sim.mat
  results$gr <- gr
  print(table(results$clustering))
  return(results)
}

cocl <- function(dm, k, f=0.66, B=3, plot.ev=TRUE, seed=42) {
  set.seed(seed)
  n <- nrow(dm)
  cl.lst <- list()
  cooc.sum <- matrix(as.integer(0), n, n)
  cosel.cnts <- matrix(as.integer(0), n, n)
  ss <- floor(n*f)
  #for (b in 1:B) {
  b <- 0
  while (min(cosel.cnts) < B) {
    b <- b + 1
    print(b)
    #cl.lst[[b]] <- rep(NA, n)
    #sel <- sample(n, floor(n*f))
    
    p <- exp(-(cosel.cnts+1))
    vals <- sample.int(length(cosel.cnts), ss*2, replace=TRUE, prob=p)
    sel <- sort(unique(c(vals %/% n + 1, vals %% n))[1:ss])
    cat('selected cells', length(sel), '\n')
    
    sim.mat <- jacc.sim(dm[sel, ], k)
    gr <- graph_from_adjacency_matrix(sim.mat, mode='undirected', weighted=TRUE, diag=FALSE)
    cl <- rep(NA, n)
    cl[sel] <- as.integer(membership(cluster_louvain(gr)))
    tmp <- outer(cl, cl, FUN='==')
    cosel.cnts <- cosel.cnts + matrix(as.integer(!is.na(tmp)), n)
    tmp[is.na(tmp)] <- as.integer(0)
    cooc.sum <- cooc.sum + tmp
    print(sum(cosel.cnts <= B-1))
    #cl.lst[[b]][sel] <- as.integer(membership(cluster_louvain(gr)))
    #tmp2 <- outer(tmp[[1]], tmp[[1]], FUN='==') + 0
  }
  cc.freq <- cooc.sum / cosel.cnts
  cc.freq[is.nan(cc.freq)] <- 0
  gr <- graph_from_adjacency_matrix(cc.freq, mode='undirected', weighted=TRUE, diag=FALSE)
  cl <- as.integer(membership(cluster_louvain(gr)))
  return(cl)
}




group9_cluster_CH1 <- function(cm) {
  #mc.cores <- 4
  #options(mc.cores = use.cores)
  cat('using', getOption("mc.cores"), 'cores\n')
  md <- data.frame(mols=apply(cm, 2, sum), genes=apply(cm>0, 2, sum))
  md$mols.per.gene <- md$mols / md$genes
  gene.th <- 200
  cat(sum(md$genes < 200), 'cells have less than 200 genes, they will be removed and get cluster label NA\n')
  sel.c <- md$genes >= 200
  regr.factors <- c('mols', 'mols.per.gene')
  gene.dr <- apply(cm[, sel.c]>0, 1, mean)
  gene.mean <- apply(cm[, sel.c], 1, mean)
  genes <- gene.dr > 0.005 & gene.mean > 0.01
  expr <- norm.nb.reg(cm[genes, sel.c], md[sel.c, regr.factors, drop=FALSE], pr.th = 30)
  cl <- cluster.the.data.simple(cm[, sel.c], expr, k=31)
  clustering <- rep(NA, ncol(cm))
  names(clustering) <- colnames(cm)
  clustering[sel.c] <- cl$clustering
  return(clustering)
}

group9_cluster_CH2 <- function(cm) {
  #options(mc.cores = use.cores)
  cat('using', getOption("mc.cores"), 'cores\n')
  md <- data.frame(mols=apply(cm, 2, sum), genes=apply(cm>0, 2, sum))
  md$mols.per.gene <- md$mols / md$genes
  gene.th <- 200
  cat(sum(md$genes < 200), 'cells have less than 200 genes, they will be removed and get cluster label NA\n')
  sel.c <- md$genes >= 200
  regr.factors <- c('mols', 'mols.per.gene')
  gene.dr <- apply(cm[, sel.c]>0, 1, mean)
  gene.mean <- apply(cm[, sel.c], 1, mean)
  genes <- gene.dr > 0.005 & gene.mean > 0.01
  expr <- norm.nb.reg(cm[genes, sel.c], md[sel.c, regr.factors, drop=FALSE], pr.th = 30)
  cl <- cluster.the.data.simple(cm[, sel.c], expr, k=31)
  cl$clustering<- cocl(cl$dm, 21, f=0.66, B=13, plot.ev=FALSE, seed=42)
  clustering <- rep(NA, ncol(cm))
  names(clustering) <- colnames(cm)
  clustering[sel.c] <- cl$clustering
  return(clustering)
}

var.genes.cm <- function(cm, min.mean=0.001, min.f.cells=0.001, z.th=1, do.plot=FALSE, span=0.33) {
  raw.mean <- apply(cm, 1, mean)
  raw.var <- apply(cm, 1, var)
  raw.sd <- apply(cm, 1, sd)
  f.cells.observed <- apply(cm > 0, 1, mean)
  sel <- f.cells.observed >= min.f.cells & raw.mean >= min.mean
  #sel <- f.cells.observed * ncol(cm) > 2
  x <- log10(raw.mean[sel])
  y <- log10(raw.var[sel]/raw.mean[sel])
  #y <- log10(raw.var[sel]/(raw.mean[sel]+0.001))
  #y <- log10(raw.sd[sel]/(raw.mean[sel]+0.001))
  fit <- loess(y ~ x, span=span, degree=2)
  gene.var.scaled <- scale(fit$residuals)[, 1]
  sel <- gene.var.scaled > z.th #& f.cells.observed[sel] >= min.f.cells & raw.mean[sel] >= min.mean
  vg <- names(which(sel))
  if (do.plot) {
    smoothScatter(x, y, xlab='log10(mean)', ylab='log10(sd/mean)')
    points(x, fit$fitted, col='green', pch=16)
    points(x[sel], y[sel], pch=16, col='deeppink')
    text(x[sel], y[sel], names(gene.var.scaled)[sel], cex=0.5)
  }
  return(vg)
}


group9_cluster_CH3 <- function(cm) {
  clustering <- rep(NA, ncol(cm))
  names(clustering) <- colnames(cm)
  gene.th <- 200
  genes <- apply(cm>0, 2, sum)
  cat(sum(genes < 200), 'cells have less than 200 genes, they will be removed and get cluster label NA\n')
  sel.c <- genes >= 200
  cm <- cm[, sel.c]
    
  vg <- var.genes.cm(cm, min.mean=0.01, min.f.cells=0.001, z.th=1, do.plot=FALSE, span=0.33)
  cat('selected', length(vg), 'variable genes\n')
  print(head(vg))
  cm.norm <- log10(cm + 1)
  cell.sum <- apply(cm.norm, 2, sum)
  cm.norm <- t(t(cm.norm[vg, ])/cell.sum)
  bpca <- batchpca(scale(cm.norm, center=TRUE, scale=TRUE), 50, byrow=TRUE)
  ev.red <- diff(bpca$values)/-bpca$values[-length(bpca$values)]
  sig.dim <- rev(which(ev.red > 0.15))[1]
  plot(bpca$values)
  abline(v=sig.dim+0.5, col='red')
  dr <- scale(bpca$vectors[, 1:sig.dim])
  js <- jacc.sim(dr, 31)
  clustering[sel.c] <- as.numeric(membership(cluster_louvain(graph_from_adjacency_matrix(js, mode='undirected', weighted=TRUE, diag=FALSE))))
  return(clustering)
}

group9_cluster_CH4 <- function(cm) {
  clustering <- rep(NA, ncol(cm))
  names(clustering) <- colnames(cm)
  gene.th <- 200
  genes <- apply(cm>0, 2, sum)
  cat(sum(genes < 200), 'cells have less than 200 genes, they will be removed and get cluster label NA\n')
  sel.c <- genes >= 200
  cm <- cm[, sel.c]
    
  vg <- var.genes.cm(cm, min.mean=0.01, min.f.cells=0.001, z.th=1, do.plot=FALSE, span=0.33)
  cm.norm <- log10(cm + 1)
  cell.sum <- apply(cm.norm, 2, sum)
  cm.norm <- t(t(cm.norm[vg, ])/cell.sum)
  bpca <- batchpca(scale(cm.norm, center=TRUE, scale=TRUE), 50, byrow=TRUE)
  ev.red <- diff(bpca$values)/-bpca$values[-length(bpca$values)]
  sig.dim <- rev(which(ev.red > 0.15))[1]
  plot(bpca$values)
  abline(v=sig.dim+0.5, col='red')
  dr <- scale(bpca$vectors[, 1:sig.dim])
  clustering[sel.c] <- cocl(dr, 21, f=0.66, B=13, plot.ev=FALSE, seed=42)
  return(clustering)
}

#load('ZeiselData.Rdata')  # "zeiselCounts"   "zeiselMetaData"
#cm <- zeiselCounts
#tmp <- group9_cluster_CH2(cm, mc.cores=4)
#tmp <- group9_cluster_CH3(cm)
#tmp <- group9_cluster_CH4(cm)
