get.Avg.Neighbor.Distance <-
function(PositionData){
  dist.matrix <- as.matrix(dist(PositionData[,4:6]))
  off.diag <- dist.matrix[col(dist.matrix)==row(dist.matrix)+1]
  return(mean(off.diag))
}
