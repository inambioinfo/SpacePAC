get.Local.Radius.Results <-
function(results.blank, position.matrix, protein.metrics, radius, culled.mutation.data, numsims, total.mutations, blank.protein, simMaxBalls){
  within.range.list<- sapply(results.blank$Center, get.distance, position.matrix[,2],protein.metrics$distance.matrix, simplify = FALSE)
  within.range <- lapply(within.range.list, function(x){as.numeric(names(which(x <radius)))} )
  results.blank$Within.Range <- within.range
  results.blank$Start <- unlist(lapply(results.blank$Within.Range,min))
  results.blank$End  <- unlist(lapply(results.blank$Within.Range,max))
  results.blank$Line.Length <- results.blank$End - results.blank$Start + 1
  
  
  
  #fill in mutation data. 
  results.obs<- calc.Sphere.Metrics(position.matrix, protein.metrics, radius, results.blank)
  results.obs$MutsCount <- unlist(results.obs$MutsCount)
  results.obs <- results.obs [order(results.obs$MutsCount, decreasing =TRUE),]
  colnames(results.blank) <- c("Line.Length","Center", "Start","End","Positions","MutsCount","Z.Score", "Within.Range")
  colnames(results.obs) <- c("Line.Length","Center", "Start","End","Positions","MutsCount","Z.Score", "Within.Range")
  results.sample <- replicate(numsims, simulate.Protein(culled.mutation.data, position.matrix, total.mutations, blank.protein, protein.metrics, radius, results.blank), simplify =FALSE)
  sim.results <- matrix(data = 0, nrow = numsims, ncol = 3)
  #Initialize all too null so we can reference them in return statement
  results.1ball <-NULL
  results.expanded2ball <-NULL
  results.expanded3ball <-NULL
  
  max.ball <-1 
  results.1ball <-  results.obs[which(results.obs$MutsCount == max(unlist(results.obs$MutsCount))),]
  max.count1Ball <- unlist(lapply(results.sample,function(x){max(unlist(x$MutsCount))}))
  one.ball.mean <- mean(max.count1Ball)
  one.ball.sd <- sd(max.count1Ball)
  if(is.na(one.ball.sd) || one.ball.sd ==0){ #This case is extremely rare and occurs only when numsims is small. However, we still account for it. 
    one.ball.sd <- 1 #in this case, every observation equals the mean. Thus when we normalize, we will get 0 as the z-score.
  }
  sim.results[,1]  <- (max.count1Ball - one.ball.mean)/one.ball.sd
  results.1ball$Z.Score <- (results.1ball$MutsCount - one.ball.mean)/one.ball.sd
  max.value <- max(results.1ball$Z.Score)
  
  if(simMaxBalls >= 2){
    if(dim(results.obs[which(results.obs$MutsCount!=0), ])[1] <2 ){ #if all muts on one position, no point in looking at 2 balls
      results.expanded2ball <- NULL
    }else{
    results.expanded2ball <- find.2Ball.Max(results.obs[which(results.obs$MutsCount!=0), ], type = "all")
      if(dim(results.expanded2ball)[1]!= 0) { #ball size ok
        max.count2Ball <- unlist(lapply(results.sample, find.2Ball.Max))
        two.ball.mean <- mean(max.count2Ball)
        two.ball.sd <-  sd(max.count2Ball)
        if(is.na(two.ball.sd) || two.ball.sd ==0){#This case is extremely rare and occurs only when numsims is small or number of muts is is very small. However, we still account for it. 
          two.ball.sd <- 1 #in this case, every observation equals the mean. Thus when we normalize, we will get 0 as the z-score.
        }
        sim.results[,2] <- (max.count2Ball - two.ball.mean)/two.ball.sd
        results.expanded2ball$Z.Score <- (results.expanded2ball$MutsCountTotal - two.ball.mean)/two.ball.sd
      
        if(max(results.expanded2ball$Z.Score) > max.value){
          max.value <- max(results.expanded2ball$Z.Score)
          max.ball <- 2
        }
      }else{
        results.expanded2ball <- NULL
      }
    }
    
    if(simMaxBalls ==3){
      if(!(is.null(results.expanded2ball))){#if 2ball has too many overlaps, no reason to do 3 ball
        if(dim(results.obs[which(results.obs$MutsCount!=0), ])[1] <3 ){ #if all muts on one position, no point in looking at 2 balls
          results.expanded3ball <- NULL
        }else{
          results.expanded3ball <- find.3Ball.Max(results.obs[which(results.obs$MutsCount!=0), ], type = "all")
          if(dim(results.expanded3ball)[1]!=0){#ball size ok
            max.count3Ball <- unlist(lapply(results.sample, find.3Ball.Max))
            three.ball.mean <- mean(max.count3Ball)
            three.ball.sd <- sd(max.count3Ball)
            if(is.na(three.ball.sd) || three.ball.sd ==0){#This case is extremely rare and occurs only when numsims is small. However, we still account for it. 
              three.ball.sd <- 1 #in this case, every observation equals the mean. Thus when we normalize, we will get 0 as the z-score.
            }
            sim.results[,3] <- (max.count3Ball - three.ball.mean)/three.ball.sd
            results.expanded3ball$Z.Score <- (results.expanded3ball$MutsCountTotal - three.ball.mean)/three.ball.sd
            
            if(max(results.expanded3ball$Z.Score) > max.value){
              max.value <- max(results.expanded3ball$Z.Score)
              max.ball <- 3
            }   
          }else{
            results.expanded3ball <- NULL
          }
        }
      }else{
      results.expanded3ball <- NULL
    }
  }
}
return(list(sim.results= sim.results, results.1ball = results.1ball, results.expanded2ball = results.expanded2ball, results.expanded3ball = results.expanded3ball, max.ball = max.ball, max.value = max.value))
}
