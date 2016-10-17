get.distance <-
function(start,destination, dist.matrix){
  row <- which(rownames(dist.matrix) == start)
  col <- which(colnames(dist.matrix) == destination)
  return(dist.matrix[row,col])
}
