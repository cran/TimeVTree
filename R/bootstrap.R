bootstrap <-
function(B=20, nodetree, subtrees, survtime, survstatus, x, 
                        D=4,minfail=30, alphac=2) {
  store <- list(NULL, NULL, NULL, NULL, NULL)
  # uses loglik (not score stat) to prune
  n <- length(survtime)
  bcoef <- store[[1]]
  btree <- store[[2]]
  bomega <- store[[3]]
  bootcost <- store[[4]]
  ori.boot <- store[[5]]
  
  K <- nrow(subtrees)
  alpha1 <- 1:K
  for (k in 1:(K-1)) alpha1[k] <- sqrt(subtrees[k,3]*subtrees[k+1,3])
  alpha1[K] <- subtrees[K,3]+1
  
  omega <- matrix(0, B, K)
  cost <- matrix(0, B, K)
  ori <- matrix(0, B, K)
  G <- omega[1,]
  
  #termblks <- nodetree[nodetree[,4] == 0,10]  
  # Number of cases in terminal blocks
  
  print('Starting Bootstrap')
  for (b in 1:B) {
    print(paste("Bootstrap Run #",b))
    #l <- numeric(0)
    #a0 <- 0
    # for (j in 1:length(termblks)) {
    # aa <- a0 + (1:termblks[j])
    # l <- c(l,sample(aa,replace=T) )
    # a0 <- a0 + termblks[j] }
    l <- sample(1:n,replace=T)
    
    l <- sort(l)
    ttime<-survtime[l]
    sstatus<-survstatus[l]
    xx<-x[l,, drop = F]
    
    # Fit the model with the bootstrap data
    boot.tree <- coxph.tree(ttime, sstatus, xx, D, minfail=minfail) 
    boot.nodetree <- output.coxphout(boot.tree)
    
    # Store bootstrap coefs and tree
    bcoef<-c(bcoef, boot.tree$coef) 
    btree<-c(btree, boot.nodetree) 
    
    # send original data thru boot.nodetree
    bootscore <- boot.coxfit(boot.nodetree,survtime,survstatus,x)
    print(bootscore)
    
    boot.subtrees <- prune(boot.nodetree)
    m <- nrow(boot.subtrees)
    
    subtree <- bootscore
    for (k in 1:K) {
      # find the optimally pruned (bootstrap) subtree for each alpha1
      i<-2
      while (i<=m) {
        if (alpha1[k]>=boot.subtrees[i,3]) {
          subtree <- free(subtree, boot.subtrees[i,5])
          i<-i+1
        }
        else break
      }
      i <- i-1
      # calculate o<-k = G(X;X<-b)-G(X<-b;X<-b) for this subtree
      G.Xb <- -boot.subtrees[i,4]
      G.X <- 0
      for (i in 1:nrow(subtree)) 
        if (subtree[i,4]!=0) {
          if (subtree[subtree[i,4], 4]==0) G.X <- G.X + subtree[subtree[i,4],7]
          if (subtree[subtree[i,5], 4]==0) G.X <- G.X + subtree[subtree[i,5],7]
        }
      if (G.X==0)   G.X<-subtree[1,7]    # the subtree is just the root
      omega[b, k] <- (G.X-G.Xb)
      cost[b, k] <- (G.Xb)
      ori[b, k] <- (G.X)
    }
    print(omega[b,])
    #cat(omega[b,], file="breast/omega", append=T)
    
  } # End of bootstrap
  
  print(omega)
  bootcost <- rbind(bootcost, cost)
  ori.boot <- rbind(ori.boot, ori)
  bomega <- rbind(bomega,omega)
  
  # find the final tree
  
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
  print(G)
  print(final)
  list( 'bcoef' = bcoef,'btree'= btree,'bomega'=  bomega,'bootcost'= bootcost,'ori.boot'= ori.boot )
}
