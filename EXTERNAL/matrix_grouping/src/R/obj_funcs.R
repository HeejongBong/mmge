obj_L0 <- function(p, L, covmat_inverse, S_bar, nn, graph, lambda1, lambda2, tau) {
  log.liklyhood <- 0
  penalty <- 0
  for (l in 1:L) {
    tmp_covmat <- covmat_inverse[,((l-1)*p+1):(l*p)]
    tmp_S <- S_bar[,((l-1)*p+1):(l*p)]
    log.liklyhood <- log.liklyhood + nn[l] * (- log(det(tmp_covmat)) 
                                + sum(tmp_covmat*tmp_S))
    tmp_covmat[abs(tmp_covmat) > tau] <- tau
    penalty <- penalty + lambda1 * (sum(apply(abs(tmp_covmat),2,sum)) - sum(abs(diag(tmp_covmat))))
  }
  graph <- graph + 1 # switch to R indexing #
  Num.of.edges <- dim(graph)[2]
  for (e in 1:Num.of.edges) {
	  l1 <- graph[1,e]
	  l2 <- graph[2,e]
	  tmp_covmat <- covmat_inverse[,(l1*p+1):((l1+1)*p)] - covmat_inverse[,((l2-1)*p+1):(l2*p)]
      tmp_covmat[abs(tmp_covmat) > tau] <- tau
      penalty <- penalty + lambda2 * (sum(apply(abs(tmp_covmat),1,sum)) - sum(abs(diag(tmp_covmat))))
      
  }
  return (log.liklyhood + penalty / 2)
}