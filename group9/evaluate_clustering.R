library("mclust")
source("~/shared_scratch/group9/task5/christoph/nmi.R")
evaluate_clustering <- function(clusts1, clusts2){
  c1 <- clusts1[intersect(names(clusts1), names(clusts2))]
  c1[is.na(c1)] <- sample(setdiff(1:(length(unique(c1))+1), unique((c1))), 1)
  c2 <- clusts2[intersect(names(clusts1), names(clusts2))]
  c2[is.na(c2)] <- sample(setdiff(1:(length(unique(c2))+1), unique((c2))), 1)
  return(list(nmi = nmi(c1, c2), ari = adjustedRandIndex(c1, c2)))
}
