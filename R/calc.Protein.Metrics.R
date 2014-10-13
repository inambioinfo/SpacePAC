calc.Protein.Metrics <-
function(culled.mutation.matrix, position.matrix){
  distance.matrix <- create.Distance.Matrix(position.matrix)
  
  #longest distance, mutated positions, and indices
  scale.factor <- max(distance.matrix)
  culled.mut.positions <- get.Mutated.Positions(culled.mutation.matrix)
  culled.mut.indices <- sapply(culled.mut.positions$names, substring, first = 2)
  
  return (list(scale.factor = scale.factor,culled.mut.counts=culled.mut.positions$counts, culled.mut.indices = culled.mut.indices, distance.matrix=distance.matrix))
}
