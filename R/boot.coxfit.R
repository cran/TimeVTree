boot.coxfit <-
function(nodetree,ttime,sstatus,x, 
                        method="breslow",iter.max=20,eps=0.000001) {
  
  # Order by survival ttime
  X <- x  ##Covariate matrix
  p <- ncol(X)  ##number of covariates
  xnames <- xname(X)
  
  n <- length(ttime)
  survindx <- 1:n
  failttime <- unique(ttime[sstatus == 1]) # unique event ttimes
  failindx <- 1:length(failttime)
  
  cutttimes <- lkl <- numeric(0)
  nevent <- nblocks <- integer(0)
  # represents number of events and blocks resp. at depth d
  bnodetree <- nodetree 
  
  bb <- 0
  
  D <- max(nodetree[,1])
  
  for (d in 0:D) {
    
    # redefine the covariates to estimate the time-varying Cox PH model
    
    cutttimesd <-  nodetree[nodetree[,1] == d,9]
    oldcutttimes <- nodetree[nodetree[,1] < d & nodetree[,4] == 0,9]
    cutttimesd <- sort(c(cutttimesd,oldcutttimes))
    
    nblocks <- c(nblocks,length(cutttimesd))
    
    breaksd <- integer(nblocks[d+1])
    for (b in 1:nblocks[d+1] ) {
      if (b == nblocks[d+1]) cutttimesd[b] <- max(ttime)
      breaksd[b] <- max(survindx[ttime <= cutttimesd[b]]) }
    breaksd <- c(0,breaksd)
    nd <- diff(breaksd) #length of blocks at depth d
    
    for (b in nodetree[nodetree[,1]==d,2]) {
      bb <- bb + 1
      blkd <- (breaksd[b]+1):breaksd[b+1]
      nnd <- (breaksd[b]+1):n
      Xdb <- X[nnd,, drop = F]
      survevent <- rep(0,length(nnd))
      survblock <-  ttime[nnd]
      survevent[1:nd[b]] <- sstatus[blkd] 
      
      if(p != 1){
        coxout <- coxph(Surv(survblock,survevent) ~ .,
                        method=method,eps=eps,iter.max=iter.max,x=T, data= Xdb)
      }else{
        run <- as.formula(paste('Surv(survblock,survevent) ~', xnames))
        coxout <- coxph(run,
                        method=method,eps=eps,iter.max=iter.max,x=T, data= Xdb)
      }
      
      names(coxout$coef) <- paste(xnames,sep=" : ",d,b)
      #print(coxout) 
      
      lkl <- c(lkl,coxout$loglik[2])
      
      bnodetree[bb,7] <- coxout$loglik[2]
      
      if (b == 1) startbl <- 0
      else startbl <- ttime[breaksd[b] + 1]
      endbl <- ttime[breaksd[b + 1]]
      
      bnodetree[bb,8] <- startbl
      bnodetree[bb,9] <- endbl
      bnodetree[bb,10] <- nd[b]
      bnodetree[bb,11] <- length(survevent[survevent==1])
    }}
  #print(lkl)
  bnodetree
}
