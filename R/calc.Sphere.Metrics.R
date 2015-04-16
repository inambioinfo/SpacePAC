calc.Sphere.Metrics <-
function(position.matrix, protein.metrics, radius, results, expanded = FALSE){

  if(!expanded){
    mut.indices.in.range<- lapply(results$Within.Range, function(x){as.numeric(protein.metrics$culled.mut.indices) %in% x})
    results$Positions <- lapply(mut.indices.in.range, function(x){as.numeric(protein.metrics$culled.mut.indices[x])}) 
    results$MutsCount <- lapply(mut.indices.in.range, function(x){ sum(protein.metrics$culled.mut.count[names(protein.metrics$culled.mut.indices[x])])})
  }else{
    mut.indices.in.range1<- lapply(results$Within.Range1, function(x){as.numeric(protein.metrics$culled.mut.indices) %in% unlist(x)})
    mut.indices.in.range2<- lapply(results$Within.Range2, function(x){as.numeric(protein.metrics$culled.mut.indices) %in% unlist(x)})
    results$Positions1 <- lapply(mut.indices.in.range1,function(x){as.numeric(protein.metrics$culled.mut.indices[x])})
    results$Positions2 <- lapply(mut.indices.in.range2,function(x){as.numeric(protein.metrics$culled.mut.indices[x])})
    results$MutsCount1 <- lapply(mut.indices.in.range1,function(x){ sum(protein.metrics$culled.mut.count[names(protein.metrics$culled.mut.indices[x])])})
    results$MutsCount2 <- lapply(mut.indices.in.range2,function(x){ sum(protein.metrics$culled.mut.count[names(protein.metrics$culled.mut.indices[x])])})
    results$MutsCountTotal <- unlist(results$MutsCount1)+unlist(results$MutsCount2)
  }
  
  return(results)
  
}
