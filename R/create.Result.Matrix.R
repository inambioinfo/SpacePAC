create.Result.Matrix <-
function(results.obs, results.sample, multcomp, alpha, method){
  
  if(method == "SimulationConservative"){
    count.results<- Reduce("+", lapply(results.sample,function(x){unlist(results.obs$MutsCount) > unlist(x$MutsCount)} ))
    p.value.column <- 1- count.results/length(results.sample)
    results.obs$P.Value<- unlist(lapply(p.value.column , function(x){ if(x == 0) {x= 1/(2*length(results.sample))} else{x = x} }))
    results.obs$P.Value <- p.adjust(results.obs$P.Value, method = multcomp)
  }else if(method == "Poisson"){    
    results.obs$P.Value <- p.adjust(results.obs$P.Value, method = multcomp)
  }
  
  sig.indices <- which(results.obs$P.Value <= alpha)
  results.obs <- results.obs[sig.indices,]
  results.obs <- results.obs[order(results.obs$P.Value,decreasing=FALSE),]
  
  return(results.obs)

}
