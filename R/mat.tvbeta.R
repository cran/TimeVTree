mat.tvbeta <-
function(indx, fulltree, subtrees = NULL, survtime, survstatus, x)
{
  Depth <- fulltree[, 1]
  Block <- fulltree[, 2]
  Node <- fulltree[, 3]
  Left <- fulltree[, 4]
  Right <- fulltree[, 5]
  Score <- -fulltree[, 7]
  Start <- fulltree[, 8]
  End <- fulltree[, 9]
  p <- ncol(x)
  xnames <- xname(x)
  #openlook()
  #postscript(file="tree1.ps")
  
  if(is.null(subtrees) == F) {
    k <- nrow(subtrees)
    subtree <- fulltree
    if (indx < k) { 
      if (indx > 1) {
        for(j in 2:indx) {
          termj <- subtrees[j, 5]
          subtree[subtree[, 3] == termj, 4] <- 0
          subtree[subtree[, 3] == termj, 5] <- 0
          nodelist <- 1:nrow(subtree)
          while(length(termj) != 0) {
            curj <- termj[length(termj)]
            termj <- termj[ - length(termj)]
            if(Left[Node == curj] != 0) {
              ll <- Left[Node == curj]
              rr <- Right[Node == curj]
              subtree <- subtree[subtree[, 3] != ll,  ]
              subtree <- subtree[subtree[, 3] != rr,  ]
              termj <- c(termj, ll, rr)
            }
          } } }
      termtree <- subtree[subtree[,4] == 0, ]
      termtree <- termtree[order(termtree[,9]),]
      #print(termtree)
      cuttimes <- sort(termtree[,9]) 
      nblocks <- length(cuttimes)
      survindx <- 1:length(survtime)
      breaks <- integer(nblocks)
      for (b in 1:nblocks) {
        breaks[b] <- max(survindx[survtime <= cuttimes[b]]) }
      breaks <- c(0,breaks)
      nd <- diff(breaks) #length of blocks 
      
      x <- as.matrix(x)
      coef <- double(ncol(x))
      
      for (b in 1:nblocks) {
        blkd <- (breaks[b]+1):breaks[b+1]
        nnd <- (breaks[b]+1):length(survtime)
        #Define start,end of block and index of failures in block
        if (b == 1) startbl <- 0
        else startbl <- survtime[breaks[b] + 1]
        endbl <- survtime[breaks[b + 1]]
        Xb <- x[nnd,, drop = F]
        Xb <- data.frame(Xb)
        survevent <- rep(0,length(nnd))
        survblock <-  survtime[nnd]
        survevent[1:nd[b]] <- survstatus[blkd]
        
        if(p != 1){
          coxout <- coxph(Surv(survblock,survevent) ~ .,
                          data= Xb)
        }else{
          run <- as.formula(paste('Surv(survblock,survevent) ~', xnames))
          coxout <- coxph(run,
                          data= Xb)
        }
        coef <- rbind(coef,coxout$coef)
        names(coxout$coef) <- paste(xnames,sep=" : ",termtree[b,3])
        print(coxout)
        
      }
      coef <- coef[-1,] 
      
      
      
      coefrep <- matrix(0,length(survtime),ncol(x))
      for (j in 1:ncol(x)) {
        coefj <- rep(coef[,j],nd)
        coefrep[,j] <- coefj  } }
    
    else if (indx == k) {
      if(p != 1){
        coxout <- coxph(Surv(survtime,survstatus) ~ .,
                        data= data.frame(x))
      }else{
        run <- as.formula(paste('Surv(survtime,survstatus) ~', xnames))
        coxout <- coxph(run,
                        data= data.frame(x))
      }
      print(coxout)
      coef <- coxout$coef
      coefrep <- matrix(0,length(survtime),ncol(x))
      for (j in 1:ncol(x)) {
        coefj <- rep(coef[j],nrow(x))
        coefrep[,j] <- coefj  } }
  }
  coefrep 
}
