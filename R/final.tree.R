final.tree <-
function(nodetree=nodetree, subtrees=subtrees, omega, alphac=2) {
  K <- nrow(subtrees)
  B <- nrow(omega)
  G <- 1:K
  
  G[1] <- -subtrees[K,4] + sum(omega[,1])/B
  G.alphac <- G[1]-alphac*subtrees[1,2]
  final <- 1
  subtree <- nodetree
  for (k in 2:K) {
    G[k] <- -subtrees[k,4] + sum(omega[,k])/B
    # next 2 lines count the number of terminal nodes of the k-th subtree
    # subtree <- free(subtree, subtrees[k,5])
    term.nodes <- subtrees[k, 2]
    if (G.alphac <= G[k]-alphac*term.nodes) {
      G.alphac <- G[k]-alphac*term.nodes
      final <- k
    }
  }
  print(G)   # this is log lik which is -cost
  print(final)
  return(list( 'subtree'= cbind(subtrees, 'cost' = (-1*(G))), 'final' = final))
}
