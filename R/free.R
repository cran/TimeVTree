free <-
function(giventree, node) {
  if (giventree[node, 4]!=0) giventree <- free(giventree, giventree[node, 4])
  if (giventree[node, 5]!=0) giventree <- free(giventree, giventree[node, 5])
  giventree[node, 4] <- 0
  giventree[node, 5] <- 0
  giventree
}
