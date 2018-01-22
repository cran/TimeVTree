prune <-
function(fulltree) {
  
  # this uses the log partial lik to prune
  
  alpha<-0
  infty<-10000
  epsilon<-.00001
  
  m<-length(fulltree[,3])
  l<-fulltree[,4]
  r<-fulltree[,5]
  #score<-fulltree[,6]
  R <- -fulltree[,7]
  N<-S<-G<-p<-g<-1:m
  subtrees<-matrix(0, m, 5)
  
  k<-1
  #for (i in m:1) {
  #   R[i] <- fulltree[i,7]
  #   if (l[i]!=0) 
  #      R[i]<-score[i]+R[l[i]]+R[r[i]]
  #}
  #R_fulltree[,7]
  
  for (i in m:1) {
    if (l[i]==0) {
      N[i]<-1
      S[i]<-R[i]
      G[i]<-infty
    }
    else {
      p[l[i]]<-p[r[i]]<-i        # memorize the parent
      N[i]<-N[l[i]]+N[r[i]]
      S[i]<-S[l[i]]+S[r[i]]
      g[i]<-(R[i]-S[i])/(N[i]-1)
      G[i]<-min(g[i], G[l[i]], G[r[i]])
    }
  }
  
  repeat {
    
    if (G[1]>alpha+epsilon) {
      subtrees[k, 1:4]<-c(k, N[1], alpha, S[1])
      alpha<-G[1]
      k<-k+1
    }
    
    if (N[1]==1)   break
    
    i<-1          #search for the weakest link
    while (G[i]<g[i]-epsilon) {
      if (G[i]==G[l[i]])   i<-l[i]
      else   i<-r[i]
    }
    subtrees[k, 5]<-i
    
    N[i]<-1
    S[i]<-R[i]
    G[i]<-infty
    while (i>1) {
      i<-p[i]     # recalculate for the tree with the branch pruned off
      N[i]<-N[l[i]]+N[r[i]]
      S[i]<-S[l[i]]+S[r[i]]
      g[i]<-(R[i]-S[i])/(N[i]-1)
      G[i]<-min(g[i], G[l[i]], G[r[i]])
    }
  }
  
  dimnames(subtrees) <- list(NULL, c('k', 'N[1]', 'alpha', 'S[1]', 'pruneoff'))
  subtrees <- subtrees[1:k-1,]
  print(subtrees)
  subtrees
}
