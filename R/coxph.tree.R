coxph.tree <- function(survtime,survstatus,x, 
                       D=3, method="breslow",minfail=10, iter.max=20,eps=0.000001, type = 'mod') {
  
  X <- x  ##Covariate matrix
  p <- ncol(X)  ##number of covariates
  n <- length(survtime)  ## number of events
  xnames <- xname(X)
  survindx <- 1:n
  failtime <- unique(survtime[survstatus == 1]) # unique event times
  failindx <- 1:length(failtime)
  type <- type
  
  beta <- lkl <- scoretest <- numeric(0)
  ntree <- integer(0)
  bb <- breakpt <- 0
  nblocks <- integer(D+1)   # represents number of blocks at depth d
  nblocks[1]  <- nodes <- 1
  nevent <- length(failtime)  # represents number of events in each block
  nodetree <- integer(10)
  mm <- floor(minfail/2)
  
  for (d in 0:D) {
    # Number of blocks prior to depth d
    if (d > 0) {bd <- sum(nblocks[1:d])   
    } else  {
      bd <- 0
    }
    
    # Figure out break points, and length of blocks at depth d
    
    breaksd <-  breakpt[bd + (1:nblocks[d+1])]
    breaksd <- c(breaksd,n)
    nodesd <- nodes[bd + (1:nblocks[d+1])]
    if (sum(nodesd) == 0) break
    nd <- diff(breaksd)         #length of blocks at depth d
    ntree <- c(ntree,nd)
    
    for (b in 1:nblocks[d+1] ) {
      blkd <- (breaksd[b]+1):breaksd[b+1]
      nnd <- (breaksd[b]+1):n   
      #Define start,end of block and index of failures in block
      if (b == 1) {startbl <- 0
      }else {
        startbl <- survtime[breaksd[b] + 1]
      }
      endbl <- survtime[breaksd[b + 1]]
      Xdb <- X[nnd, , drop = F]
      survevent <- rep(0,length(nnd))
      survblock <- survtime[nnd]
      survevent[1:nd[b]] <- survstatus[blkd]  
      
      if(p != 1){
        coxout <- coxph(Surv(survblock,survevent) ~ .,
                        method=method,eps=eps,iter.max=iter.max,x=T, data= Xdb)
      }else{
        run <- as.formula(paste('Surv(survblock,survevent) ~', xnames))
        coxout <- coxph(run,
                        method=method,eps=eps,iter.max=iter.max,x=T, data= Xdb)
      }
      
      
      names(coxout$coef) <- paste(names(coxout$coef),sep=" : ",d,b) 
      
      print(coxout)
      
      beta <- c(beta,coxout$coef)
      coxoutd <- coxph.detail(coxout)
      scoreresid <- apply(coxoutd$score,2,rawscore) 
      if(type == 'mod'){
        infscore <- apply(coxoutd$imat,c(1,2),sum)
        infsvd <- svd(infscore)
        infpre <- infsvd$u %*% diag(1/sqrt(infsvd$d))
        #infinv <- solve(infscore)
        
        failblock <- unique(survblock[survevent == 1])
        failind <- 1:length(failblock)
        
        # Determine breakpoints
        
        scoremat <- apply(scoreresid,1,'infmul',infpre)
        scorevec <- apply(scoremat^2,2,sum)
        scoreval <- max(scorevec)
        
        scoretest <- c(scoretest,scoreval)
      }else if(type == 'ori'){
        inftot <- apply(coxoutd$imat,c(1,2),sum)
        infinv <- solve(inftot)
        infpre <- apply(coxoutd$imat,c(1,2),cumsum)
        
        failblock <- unique(survblock[survevent == 1])
        failind <- 1:length(failblock)
        
        scorevec <- rep(0,length(failblock))
        # Determine breakpoints
        
        for (jj in mm:(length(failblock)-mm)) {
          scoremat <- infpre[jj, ,] - infpre[jj, ,]%*%infinv%*%infpre[jj, ,]
          scoreinf <- solve(scoremat)
          scorevec[jj] <- t(scoreresid[jj,])%*%scoreinf%*%scoreresid[jj,] 
        }
        scoreval <- max(scorevec)
        scoretest  <- c(scoretest,scoreval)
      } else{
        warning('Unknown type was entered into the argument.')
      }
      tt <- order(scorevec)
      lftbreakpt <- breakpt[bd+b]
      rgtbreakpt <- failblock[tt[length(tt)]] 
      fleftblk <- failind[failblock <= rgtbreakpt]
      frightblk <- failind[failblock > rgtbreakpt]
      # No. of events in children of block(b,d)
      neventleft <- length(fleftblk)
      neventrght <- length(frightblk)
      
      if (((neventleft >= minfail & neventrght >= minfail) | 
           (scoreval < .0001)) & (nodes[bd+b] == 1)) {
        lkl <- c(lkl,coxout$loglik[2])
        nevent <- c(nevent,neventleft,neventrght) 
        rightbreak <- max((1:n)[survtime <= rgtbreakpt])
        breakpt <- c(breakpt,lftbreakpt,rightbreak) 
        nblocks[d+2] <- nblocks[d+2] + 2
        nodes <- c(nodes,1,1) 
        bb <- bb + 1
        if (d == 0) { 
          leftchild <- 2 
        } else {
          leftchild <- max(nodetree[,5]) + 1 
        }
        rightchild <- leftchild + 1
        if  (d == D) {
          leftchild <- rightchild <- 0
        }
        nodetree <- rbind(nodetree,c(d,b,bb,
                                     leftchild,rightchild,scoreval,
                                     #breaksd[b],breaksd[b+1],
                                     startbl,endbl,
                                     nd[b],nevent[bd+b])) 
      } else { 
        breakpt <- c(breakpt,lftbreakpt) 
        print(neventleft)
        print(neventrght)
        print(scoreval)
        nblocks[d+2] <- nblocks[d+2] + 1 
        nevent <- c(nevent,neventleft+neventrght) 
        nodes <- c(nodes,0)
        if  (nodes[bd+b] == 1) { 
          lkl <- c(lkl,coxout$loglik[2])
          bb <- bb + 1
          nodetree <- rbind(nodetree,c(d,b,bb,
                                       0,0,scoreval,
                                       #breaksd[b],breaksd[b+1],
                                       startbl,endbl,
                                       nd[b],nevent[bd+b]))} }
    }
  }
  breakpt[breakpt!= 0] <- survtime[breakpt[breakpt!= 0]]
  nodetree[nodetree[,2] == 1,7]  <-  0
  coxphtree <- list(D=D,coef=beta,lkl=lkl,breakpt=breakpt,
                    ntree=ntree,nevent=nevent,nblocks=nblocks,
                    nodes=nodes,nodetree=nodetree,scoretest=scoretest,
                    xnames=xnames,
                    #tvarind=tvarind,
                    failtime=failtime) 
  coxphtree
}
