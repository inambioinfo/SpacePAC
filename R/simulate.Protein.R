simulate.Protein <-
function(culled.mutation.data, position.matrix, total.mutations, blank.protein, protein.metrics, radius, results, expanded = FALSE){
  #Generate a random sample
  random.sample <- table(sample(x=dim(culled.mutation.data)[2], size = total.mutations, replace= TRUE ))
  simulated.protein <- blank.protein
  simulated.protein[as.numeric(names(random.sample))]<- random.sample
  
  #adjusts the relevant mutational distribution stored in protein metrics
  sim.protein.metrics <- protein.metrics
  sim.protein.metrics$culled.mut.counts <- simulated.protein
  sim.protein.metrics$culled.mut.indices <- sapply(names(which(simulated.protein > 0)), substring, first = 2)
  
  #creates the results matrix
  result <- calc.Sphere.Metrics(position.matrix, sim.protein.metrics, radius, results, expanded)
  return (result)
  
}
