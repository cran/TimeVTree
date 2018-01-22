output.coxphout <-
function(coxout) { 
  
  nodetree <- coxout$nodetree[-1,]
  Lkl <- coxout$lkl
  Depth <- nodetree[,1]
  Block <- nodetree[,2]
  Node <- nodetree[,3]
  Left <- nodetree[,4]
  Right <- nodetree[,5]
  Score <- nodetree[,6]
  Start <- nodetree[,7]
  End <- nodetree[,8]
  ncases <- nodetree[,9]
  nevents <- nodetree[,10]
  nnodetree <- as.data.frame(cbind(Depth,Block,Node,Left,
                                   Right,Score,Lkl,Start,End,ncases,nevents))
  dimnames(nnodetree) <- list(nodetree[,3],c('Depth','Block','Node','Left',
                                             'Right','Score','lkl','Start','End','# of Cases','# of Events'))
  #print(nnodetree)
  #print(coxout$coef)
  nnodetree }
