library("osDesign")
library("plyr")
library("doMC")
registerDoMC(4)

ambient_noise <- function(umis, epsilon = .05){
  n_genes <- nrow(umis)
  n_cells <- ncol(umis)
  migrating_umis <- aaply(1:nrow(umis), 1, function(g){
    rmvhyper(as.vector(umis[g,]), rpois(1, sum(umis[g,]) * epsilon))
  }, .parallel = F)
  returning_umis <- t(sapply(rowSums(migrating_umis), function(k) rmultinom(1, k, rep(1, n_cells))))
  res <- umis - migrating_umis + returning_umis
  rownames(res) <- rownames(umis)
  colnames(res) <- colnames(umis)
  return(res)
}

subsample <- function(umis, rate = .8){
  res <- Matrix(t(aaply(1:ncol(umis), 1, function(cell){
    return(rmvhyper(as.vector(umis[, cell]), rpois(1, sum(umis[, cell]) * rate)))
  }, .parallel = T)), sparse = T)
  rownames(res) <- rownames(umis)
  colnames(res) <- colnames(umis)
  return(res)
}

sample_cells <- function(umis, rate){
  return(umis[, sample(1:ncol(umis), round(ncol(umis) * rate))])
}

