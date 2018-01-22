elbow.tree <-
function(nodetree=nodetree, subtrees=subtrees, omega, alphac=2) {
  K <- nrow(subtrees)
  B <- nrow(omega)
  G <- 1:K
  G.p <- 1:K
  
  G[1] <- -subtrees[K,4] + sum(omega[,1])/B
  G.p[1] <- G[1]-alphac*subtrees[1,2]
  subtree <- nodetree
  for (k in 2:K) {
    G[k] <- -subtrees[k,4] + sum(omega[,k])/B
    # next 2 lines count the number of terminal nodes of the k-th subtree
    # subtree <- free(subtree, subtrees[k,5])
    term.nodes <- subtrees[k, 2]
    G.p[k] <- G[k]-alphac*term.nodes
  }
  return(list( 'subtree'= cbind(subtrees, 'cost' = (-1*(G)), 'cost.p' = (-1*(G.p)))))
}
