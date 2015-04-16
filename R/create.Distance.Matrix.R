create.Distance.Matrix <-
function(position.matrix){
  distance.matrix <- as.matrix(dist(position.matrix[,4:6], diag = TRUE, upper = FALSE))
  colnames(distance.matrix) <- position.matrix[,2]
  rownames(distance.matrix) <- position.matrix[,2]
  return(distance.matrix)
  
}
