get.Mutated.Positions <-
function (mutation.matrix){
  return (list(names = names(which(colSums(mutation.matrix) >0))
               , counts =  colSums(mutation.matrix)))
  
}
