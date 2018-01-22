plot_coxphtree <-
  function(fulltree, subtrees = NULL, mm = 3, start = 0, pdf = FALSE, file.name) {
    Depth <- fulltree[, 1]
    Block <- fulltree[, 2]
    Node <- fulltree[, 3]
    Left <- fulltree[, 4]
    Right <- fulltree[, 5]
    Score <- -fulltree[, 7]
    Start <- fulltree[, 8]
    End <- fulltree[,9] 
    #Score _ (Score - max(Score))/(min(Score)-max(Score))
    if(is.null(dev.list()) == F)
      graphics.off()
    #openlook()
    #ps.options(colormodel = "gray") 
    #postscript(file="tree0.ps")
    for(i in 1:length(Start)){
      if(Start[i] == 0){
        Start[i] <- start
      }
    }
    if(pdf == TRUE){
      pdf(paste(file.name, '.pdf', sep = ''))
    }
    Dmax <- 2 * max(Depth)	# Plotting the fulltree 
    # Create empty plot
    par(col=5,lwd=3)
    plot(c(start, End[1]), c(0, Dmax + 1), type = "n", axes = F, xlab = 
           "Survival Time", ylab = " ")
    axis(side = 1, at = start:floor(max(End)))	#Plot Tree
    for(i in 1:nrow(fulltree)) {
      if(Left[i] == 0) coli <- 6
      else coli <- par("col")
      
      #Terminal Node 
      polygon(c(Start[i], Start[i], End[i], End[i]), c(Dmax - 2 * 
                                                         Depth[i], Dmax - 2 * Depth[i] + 1, Dmax - 2 * Depth[i] + 
                                                         1, Dmax - 2 * Depth[i]), density = 0, col = coli)
      if (Left[i] == 0) {
        text(0.5 * (End[i] + Start[i]), Dmax - 2 * Depth[i] +
               0.5, paste(i),col=coli) 
      }
      else {text(0.5 * (End[i] + Start[i]), Dmax - 2 * Depth[i]
                 + 0.5,paste(i))}
    }
    #openlook()	# Calculate Tree heights and put in yd as you go along
    yd <- rep(0,nrow(fulltree))
    for(i in 1:nrow(fulltree)) {
      if(Left[i] != 0) {
        yd[Left[i]] <- yd[i] - Score[i] + (Score[Left[i]] + 
                                             Score[Right[i]])
        yd[Right[i]] <- yd[Left[i]]
      }
    }
    
    par(col=5,lwd=3)
    plot(c(0, 1), c(0, min(yd)), type = "n", axes = F, xlab = 
           "Height of Branch = Magnitude of Change in
         Log-Likelihood", ylab = " ")	
    #Plot Tree
    noderows <- 1:nrow(fulltree)
    for(d in 0:max(Depth)) {
      nodecurd <- noderows[(Depth < d & Left == 0) | (Depth == d)]
      nodeord <- order(Start[nodecurd])
      nodecurd <- nodecurd[nodeord]
      bb <- 0
      for(b in nodecurd) {
        bb <- bb + 2^(d - Depth[b])
        wd <- 0.5/(2^d)
        xd <- wd + 2 * (bb - 1) * wd
        if(Left[b] != 0) {
          text(xd + 0.01, yd[b] + 0.2, paste(b))
          lines(c(xd - wd/2, xd + wd/2), c(yd[b], yd[b]), col = 4)
          lines(c(xd - wd/2, xd - wd/2), c(yd[b], yd[Left[b]]), col = 4)
          lines(c(xd + wd/2, xd + wd/2), c(yd[b], yd[Right[b]]), col = 4)
        }
        else if(Depth[b] == d) {
          text(xd, yd[b] - 0.2, paste(b))
        }
      }
    }
    
    if(is.null(subtrees) == F) {
      plotno <- 1
      #openlook()
      par(mfrow = c(mm, 2),col=5,lwd=2)
      Dmax <- 2 * max(Depth)
      k <- nrow(subtrees)
      subtree <- fulltree
      for(j in 2:(k - 1)) {
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
        }
        subDepth <- subtree[, 1]
        subBlock <- subtree[, 2]
        subNode <- subtree[, 3]
        subLeft <- subtree[, 4]
        subRight <- subtree[, 5]
        subScore <- -subtree[, 7]
        subStart <- subtree[, 8]
        subEnd <- subtree[, 9]	
        for(i in 1:length(subStart)){
          if(subStart[i] == 0){
            subStart[i] <- start
          }
        }
        # Calculate Tree heights and put in yd as you go along
        subyd <- rep(0, nrow(subtree))
        for(i in 1:nrow(subtree)) {
          if(subLeft[i] != 0) {
            subyd[subNode == subLeft[i]] <- subyd[i] - subScore[i]
            subyd[subNode == subRight[i]] <- subyd[i] - subScore[i]
          }
        }
        subyd <- (subyd/subyd[[2]])* -10
        # Plotting the fulltree 
        # Create empty plot
        plot(c(start, End[1]), c(0, Dmax + 1), type = "n", axes = F,
             xlab = "Survival Time", ylab = " ")
        axis(side = 1, at = start:floor(max(End)))	#Plot Tree
        for(i in 1:nrow(subtree)) {
          if(subLeft[i] == 0) 
            coli <- 6 else coli <- par("col")
            polygon(c(subStart[i], subStart[i], subEnd[i], 
                      subEnd[i]), c(Dmax - 2 * subDepth[i], Dmax - 
                                      2 * subDepth[i] + 1, Dmax - 2 * subDepth[i] + 
                                      1, Dmax - 2 * subDepth[i]), density = 0, col
                    = coli)
            text(0.5 * (subEnd[i] + subStart[i]), Dmax - 2 * 
                   subDepth[i] + 0.5, paste(subNode[i]), cex = 
                   0.5)
        }
        
        plot(c(0, 1), c(0, min(yd)), type = "n", axes = F, xlab = " ", ylab = " ")	#Plot Tree
        title(paste('subtree', j, sep = ' '), cex.main = 2.5)
        noderows <- 1:nrow(subtree)
        for(d in 0:max(subDepth)) {
          nodecurd <- noderows[(subDepth < d & subLeft == 0) | (subDepth == d)]
          nodeord <- order(subStart[nodecurd])
          nodecurd <- subtree[nodecurd[nodeord], 3]
          bb <- 0
          for(b in nodecurd) {
            bb <- bb + 2^(d - subDepth[subNode == b])
            wd <- 0.5/(2^d)
            xd <- wd + 2 * (bb - 1) * wd	
            #print(c(d,b))
            #print(c(xd,wd))
            #print(subyd[subNode == b])
            if(subLeft[subNode == b] != 0) {
              #print(subLeft[subNode == b])
              #print(subyd[subNode == subLeft[subNode == b]])
              #print(subRight[subNode == b])
              #print(subyd[subNode == subRight[subNode == b]])
              text(xd + 0.02, subyd[subNode == b], paste(b), cex = 0.5)
              lines(c(xd - wd/2, xd + wd/2), c(subyd[subNode == b], subyd[subNode == b]), col = 4)
              lines(c(xd - wd/2, xd - wd/2), c(subyd[ subNode == b], subyd[subNode == subLeft[subNode == b]]), col = 4)
              lines(c(xd + wd/2, xd + wd/2), c(subyd[subNode == b], subyd[subNode == subRight[subNode == b]]), col = 4)
            }
            else if(subDepth[subNode == b] == d) {
              text(xd, subyd[subNode == b] - 1.5, paste(b), cex = 0.5)
            }
          }
        }
        plotno <- plotno + 1
        #dev.ask()
      }
      # Plotting the last pruned tree with 1 node
      # Create empty plot
      plot(c(start, End[1]), c(0, Dmax + 1), type = "n", axes = F, xlab
           = "Survival\nTime", ylab = " ")
      axis(side = 1, at = start:floor(max(End)))	#Plot Tree
      coli <- 6	#Terminal Node
      polygon(c(Start[1], Start[1], End[1], End[1]), c(Dmax, Dmax + 1,
                                                       Dmax + 1, Dmax), density = 0, col = coli)
      text(0.5 * (End[1] + Start[1]), Dmax + 0.5, paste(Node[1]), cex
           = 0.5)
    }
    
    if(pdf == TRUE){
      dev.off()
    }
  }
