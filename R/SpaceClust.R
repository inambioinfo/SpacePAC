SpaceClust <-
function(mutation.data, position.matrix, method = "SimMax", numsims = 1000,  simMaxSpheres = 3, radii.vector, multcomp = "bonferroni", alpha = .05){
  
  #calculates basic mutation protein information
  culled.mutation.data <- cull.Mutation.Matrix(mutation.data,position.matrix)
  protein.metrics <- calc.Protein.Metrics(culled.mutation.data, position.matrix)
  
  total.mutations <- sum(culled.mutation.data)
  blank.protein <- protein.metrics$culled.mut.counts
  blank.protein[1:length(blank.protein)] <-0
  
  #sets up the results matrix
  results.blank <- matrix(data = 0, nrow = dim(position.matrix)[1], ncol = 8)
  colnames(results.blank) <- c("Line.Length","Center", "Start","End","Positions","MutsCount","P.Value", "Within.Range")
  results.blank <- as.data.frame(results.blank)
  results.blank$Center <- position.matrix[,2]
  
  if(total.mutations > 1){
    
    if(method == "SimulationConservative"){
      results <- list()
      for(i in 1:length(radii.vector)){
        radius <- radii.vector[i]
        cat(sprintf("Processing radius # %d : radius length = %0.2f : Percentage complete %.2f \n", i, radius, (i-1)/length(radii.vector)))
        flush.console()
        within.range.list<- sapply(results.blank$Center, get.distance, position.matrix[,2],protein.metrics$distance.matrix, simplify = FALSE)
        within.range <- lapply(within.range.list, function(x){as.numeric(names(which(x <radius)))} )
        results.blank$Within.Range <- within.range
        results.blank$Start <- unlist(lapply(results.blank$Within.Range,min))
        results.blank$End  <- unlist(lapply(results.blank$Within.Range,max))
        results.blank$Line.Length <- results.blank$End - results.blank$Start + 1
      
        #fill in mutation data. Drop all balls that have 0 mutations.
        results.obs<- calc.Sphere.Metrics(position.matrix, protein.metrics, radius, results.blank)
        results.obs$MutsCount <- unlist(results.obs$MutsCount)
        
        results.sample <- replicate(numsims, simulate.Protein(culled.mutation.data, position.matrix, total.mutations, blank.protein, protein.metrics, radius, results.blank), simplify =FALSE)
        results[[i]] <- create.Result.Matrix(results.obs=results.obs, results.sample=results.sample, multcomp=multcomp,alpha=alpha, method=method)
      }
      p.val.list <- unlist(lapply(results, function(x){x$P.Value[1]})) #If two radii have the same min.pvalue, picks the smaller
      if(!all(is.na(p.val.list))){
        best.p.val <- which(p.val.list == min(p.val.list,na.rm=TRUE))[1] #If some radii are NA (meaning no clusters, they will be ignored)
        results <- results[[best.p.val]]
        results$P.Value <- results$P.Value * length(radii.vector)
        results <- list(result.matrix = results, best.radii = radii.vector[best.p.val])
      }else{
        results <- list(result.matrix = results[[1]], best.radii = NULL)
      }
    }else if(method == "Poisson"){
      results <- list()
      for(i in 1:length(radii.vector)){
        radius <- radii.vector[i]
        cat(sprintf("Processing radius # %d : radius length = %0.2f : Percentage complete %.2f \n", i, radius, (i-1)/length(radii.vector)))
        flush.console()
        within.range.list<- sapply(results.blank$Center, get.distance, position.matrix[,2],protein.metrics$distance.matrix, simplify = FALSE)
        within.range <- lapply(within.range.list, function(x){as.numeric(names(which(x <radius)))} )
        results.blank$Within.Range <- within.range
        results.blank$Start <- unlist(lapply(results.blank$Within.Range,min))
        results.blank$End  <- unlist(lapply(results.blank$Within.Range,max))
        results.blank$Line.Length <- results.blank$End - results.blank$Start + 1
        
        #fill in mutation data. Drop all balls that have 0 mutations.
        results.obs<- calc.Sphere.Metrics(position.matrix, protein.metrics, radius, results.blank)
        results.obs$MutsCount <- unlist(results.obs$MutsCount)
        lambda <- total.mutations/ dim(position.matrix)[1]
        scaled.lambdas <- lapply(results.obs$Within.Range, function(x){length(x)*lambda})
        
        #ppois finds the cumulative probability of being greater than observed count. We want greater than or equal so I subtract one.
        #IE, p-value if observed = 50 , is Pr(X >=50). PPois finds Pr(X >50) so I calculate P(X>49) = Pr(X>=50)
        p.values <- ppois(unlist(results.obs$MutsCount)-1,lambda = unlist(scaled.lambdas), lower.tail= FALSE)
        results.obs$P.Value <- p.values
        results[[i]] <- create.Result.Matrix(results.obs = results.obs, results.sample = NULL, multcomp=multcomp,alpha=alpha, method=method) 
      }
      p.val.list <- unlist(lapply(results, function(x){x$P.Value[1]})) #If two radii have the same min.pvalue, picks the smaller
      if(!all(is.na(p.val.list))){
        best.p.val <- which(p.val.list == min(p.val.list,na.rm=TRUE))[1] #If some radii are NA (meaning no clusters, they will be ignored)
        results <- results[[best.p.val]]
        results$P.Value <- results$P.Value * length(radii.vector)
        results <- list(result.poisson = results, best.radii = radii.vector[best.p.val])
      }else{
        results <- list(result.poisson = results[[1]], best.radii = NULL)
      }
    }else if(method =="SimMax"){  
      if(simMaxSpheres > 3 || simMaxSpheres < 1){
        stop("Error! The maximum number of spheres should be either 1, 2, or 3.")
      }else{
        radius.results <- list()
        best.1ball = NULL
        best.2ball = NULL
        best.3ball = NULL
        best.1ball.radius = NULL 
        best.2ball.radius = NULL
        best.3ball.radius = NULL
        bad.2ball.message = NULL
        bad.3ball.message = NULL
        bad.2ball.radii = NULL
        bad.3ball.radii = NULL
        
        for(i in 1 :length(radii.vector)){
          radius = radii.vector[i]
          cat(sprintf("Processing radius # %d : radius length = %0.2f : Percentage complete %.2f \n", i, radius, (i-1)/length(radii.vector)))
          flush.console()
          radius.results[[i]] <- get.Local.Radius.Results(results.blank, position.matrix, protein.metrics, radius, culled.mutation.data,numsims, total.mutations, blank.protein,simMaxSpheres)
          gc()
        }
        #gets the max over all radii and ball sizes
        
        #gets the max simulation result for each radius size. This provides a matrix with the number of columsn equal to the number of radii considered.
        #First the max within each readius is obtained. Then we max over all the radii.
        max.sim.results  <- lapply(radius.results, function(x){apply(x$sim.results, 1, max)})
        max.matrix <- apply(matrix(data = unlist(max.sim.results), ncol = length(max.sim.results)), 1, max) #collapses over radii
        
        #checks to see if spheres small enough
        num.balls.considered <- unlist(lapply(radius.results, function(x){ dim(x$results.1ball)[1]}))[1] #will equal total if every ball overlaps
        if(num.balls.considered == dim(position.matrix)[1]){
          stop("Error! Radii to big. Every sphere touches every other position for at least one of your radii! Reduce the size of your radii!")
        }else{
          
          #gets the scores for 1 ball, 2 ball and 3ball
          #gets which radius is best for each ball
          #gets the bestball.data
          #conditions appropriately based upon the max number of balls selected.
          z.scores.1ball <- unlist(lapply(radius.results, function(x){max(x$results.1ball$Z.Score)}))
          best.1ball.radius <- radii.vector[which(z.scores.1ball == max(z.scores.1ball))][1] #in degenerate case that several radii have same z-score, minimum radius chosen
          best.1ball <- radius.results[[which(z.scores.1ball == max(z.scores.1ball))[1]]]$results.1ball #in degenerate case that several radii have same z-score, minimum radius chosen
          max.score <- max(z.scores.1ball[which(z.scores.1ball == max(z.scores.1ball))[1]]) #in degenerate case that several radii have same z-score, minimum radius chosen
          
          #see if any balls came out unintentionally as blank.
          #this would occur if a user specifies too many balls given mutations
          #for instance, if there were mutations at just positions 10 and 11 and the user specifies 3 balls, this would not work.
          
          
          #bad radii1ball should never happen since we check if mutation count >0
          bad.radii1ball <- which(unlist(lapply(lapply(radius.results, function(x){x$results.1ball}), function(x){  if(is.null(x)) {1} else{0}  })) != 0)
          if(simMaxSpheres >= 2){
            bad.radii2ball <- which(unlist(lapply(lapply(radius.results, function(x){x$results.expanded2ball}), function(x){  if(is.null(x)) {1} else{0}  })) != 0)
            if(length(bad.radii2ball)!=0){
              bad.2ball.message = "Error when looking at 2 spheres. Most likely caused by too many balls given mutation positions. If all mutations on one position then you have this error. Can also occur if all valid 2-sphere combinations have no mutations within them. Radii at which error occured showed in 'bad.2ball.radius'."
              bad.2ball.radii <- radii.vector[bad.radii2ball]
            }
            if(simMaxSpheres ==3){
              bad.radii3ball <- which(unlist(lapply(lapply(radius.results, function(x){x$results.expanded3ball}), function(x){  if(is.null(x)) {1} else{0}  })) != 0)
              if(length(bad.radii3ball)!=0){
                bad.3ball.message = "Error when looking at 3 spheres. Most likely caused by too many balls given mutation positions. If all mutations on two position then you have this error. Can also occur if all valid 2-sphere combinations have no mutations within them. Radii at which error occured showed in 'bad.3ball.radius'."
                bad.3ball.radii <- radii.vector[bad.radii3ball]
              }
              
            }  
          }
          
          
          if(simMaxSpheres >= 2 && length(bad.radii2ball)==0){
            z.scores.2ball <- unlist(lapply(radius.results, function(x){max(x$results.expanded2ball$Z.Score)}))
            best.2ball.radius <- radii.vector[which(z.scores.2ball == max(z.scores.2ball))][1] #in degenerate case that several radii have same z-score, minimum radius chosen
            best.2ball <- radius.results[[which(z.scores.2ball == max(z.scores.2ball))[1]]]$results.expanded2ball
            max.score <- max(z.scores.1ball[which(z.scores.1ball == max(z.scores.1ball))[1]], z.scores.2ball[which(z.scores.2ball == max(z.scores.2ball))])
          }
          if(simMaxSpheres >=3 && length(bad.radii3ball)==0){
            z.scores.3ball <- unlist(lapply(radius.results, function(x){max(x$results.expanded3ball$Z.Score)}))
            best.3ball.radius <- radii.vector[which(z.scores.3ball == max(z.scores.3ball))][1]#in degenerate case that several radii have same z-score, minimum radius chosen
            best.3ball <- radius.results[[which(z.scores.3ball == max(z.scores.3ball))[1]]]$results.expanded3ball
            max.score <- max(z.scores.1ball[which(z.scores.1ball == max(z.scores.1ball))[1]], z.scores.2ball[which(z.scores.2ball == max(z.scores.2ball))], z.scores.3ball[which(z.scores.3ball == max(z.scores.3ball))])
          }
          
          if(max.score %in% z.scores.1ball){
            radius.of.max.score <- radii.vector[which(z.scores.1ball == max.score)][1]
            best.ball.count <- 1
            best.radius <- which(z.scores.1ball == max.score)[1]
            best.sphere <- radius.results[[best.radius]]$results.1ball
          }else if(simMaxSpheres >= 2 && max.score %in% z.scores.2ball && length(bad.radii2ball)<length(radii.vector) ){
            radius.of.max.score <- radii.vector[which(z.scores.2ball == max.score)][1]
            best.ball.count <-2 
            best.radius <- which(z.scores.2ball == max.score)[1]
            best.sphere <- radius.results[[best.radius]]$results.expanded2ball
          }else if(simMaxSpheres >= 3 && max.score %in% z.scores.3ball && length(bad.radii3ball)<length(radii.vector)){
            radius.of.max.score <- radii.vector[which(z.scores.3ball == max.score)][1]
            best.ball.count <-3
            best.radius <- which(z.scores.3ball == max.score)[1]
            best.sphere <- radius.results[[best.radius]]$results.expanded3ball
          }
          
          p.value <- 1 - sum(max.score > max.matrix)/numsims
          if(p.value ==0){
            p.value = 1/(2*numsims)
          }else{
            p.value = p.value
          }  
          results <- list(p.value = p.value, optimal.num.spheres = best.ball.count, optimal.radius = radius.of.max.score, optimal.sphere = best.sphere, 
                          best.1.sphere =best.1ball, best.2.sphere = best.2ball, best.3.sphere = best.3ball, best.1.sphere.radius = best.1ball.radius, 
                          best.2.sphere.radius = best.2ball.radius, best.3.sphere.radius = best.3ball.radius, bad.2.sphere.message = bad.2ball.message ,
                          bad.3.sphere.message = bad.3ball.message ,bad.2.sphere.radii = bad.2ball.radii, bad.3.sphere.radii = bad.3ball.radii)
        }
        
        
      }
    }else{
      stop("Error! Method Undetected! Must be `SimulationConservative', 'Poisson' or 'SimMax'.")
    }
  }else{
    stop("Error! Either 0 or 1 mutation in the culled data! No clusters possible!")
  }
  return(results)
}
