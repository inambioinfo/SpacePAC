find.2Ball.Max <-
function(sample, type = "first"){
  total.found = 0
  
  sample$MutsCount <- unlist(sample$MutsCount)
  sample <- sample[order(sample$MutsCount,decreasing=TRUE),]


  v<- sample$MutsCount
  l<- length(v)
  
  starti=2
  startj=1
  
  max.matrix <- matrix(data = 0, nrow = 0, ncol = 3)
  current.max <- 0
  colnames(max.matrix)<- c("i","j", "sum")
  
  cand<- as.data.frame(matrix(data = 0, nrow = 0, ncol =3))
  cand <- rbind(cand, c(v[starti] + v[startj],starti,startj), deparse.level =0)
  names(cand)<- c("sum","i","j")
  
  while(dim(cand)[1]!= 0){
    index <- which.max(cand$sum)
    sum  <- cand$sum[index]
    i <- cand$i[index]
    j <- cand$j[index]
    cand <- cand[-index,]

    intersection <- intersect(unlist(sample$Within.Range[i]), unlist(sample$Within.Range[j]))
    if ((!(i==j)) & (length(intersection) ==0)){  
      #cat(sprintf("(%d,%d,%d)\n", sum, i, j))
      total.found <- total.found +1
      if(total.found ==1){
        current.max <- sample$MutsCount[i]+sample$MutsCount[j]
      }
      
      if(type == "first"){
        break
      }else{
        if(sample$MutsCount[i]+sample$MutsCount[j] == current.max){ #still finding equivalent size balls
          max.matrix <- rbind(max.matrix, c(sample$Center[i],sample$Center[j],sample$MutsCount[i]+sample$MutsCount[j] ))
          total.found <- total.found +1
        }
        else{
          break
        }
      }
    }   
    
    if(j ==1 & i < l){
      cand <- rbind(cand, c(v[i+1]+ v[j],i+1,j), deparse.level =0)
      names(cand)<- c("sum","i","j")
    }
    
    if(j+1 < i){
      cand <- rbind(cand,  c(v[i]+ v[j+1], i, j+1),  deparse.level =0)
      names(cand)<- c("sum","i","j")
    }
  }
    if(total.found == 0){ #This is an error case, all 2 balls overlap
      max.covered = -1
    }else{ #A success was found
      max.covered = sample$MutsCount[i]+sample$MutsCount[j]
    }
  
  if(type == "first"){
    return(max.covered)
  }else{
    if(max.covered!=-1){
      results.expanded <- matrix(data = 0, nrow = dim(max.matrix)[1], ncol = 17)
      results.expanded <- as.data.frame(results.expanded)
      colnames(results.expanded) <- c("Line.Length1","Line.Length2", "Center1", "Center2", "Start1","End1","Start2", "End2", "Positions1","Positions2","MutsCount1","MutsCount2","MutsCountTotal", "Z.Score", "Within.Range1", "Within.Range2","Intersection")
      results.expanded$Center1 <- max.matrix[,1]
      results.expanded$Center2 <- max.matrix[,2]
      results.expanded$MutsCountTotal <- max.matrix[,3]
    
      part1 <- lapply(results.expanded$Center1, function(x){sample[which(x == sample$Center) ,c("Line.Length","Start","End","Positions","MutsCount","Within.Range")]})
      part2 <- lapply(results.expanded$Center2, function(x){sample[which(x == sample$Center) ,c("Line.Length","Start","End","Positions","MutsCount","Within.Range")]})
      results.expanded[,c("Line.Length1","Start1","End1", "Positions1", "MutsCount1","Within.Range1")] <-  as.data.frame(t(sapply(part1,function(x){x})))
      results.expanded[,c("Line.Length2","Start2","End2", "Positions2", "MutsCount2","Within.Range2")] <-  as.data.frame(t(sapply(part2,function(x){x})))
      
      #finds the ones with no intersection
      results.expanded$Intersection <- lapply(1:dim(results.expanded)[1], function(x){intersect(unlist(results.expanded$Within.Range1[x]), unlist(results.expanded$Within.Range2[x]))})
      no.intersect.indices <- which( unlist(lapply(results.expanded$Intersection, function(x){length(x)}))==0)
      results.expanded <- results.expanded[no.intersect.indices,]
    }
    else{#if all 2 balls overlap, return a matrix with 0 rows which will be tested for later.
      results.expanded <- matrix(data = 0, nrow = 0, ncol = 17) 
      results.expanded <- as.data.frame(results.expanded)
      colnames(results.expanded) <- c("Line.Length1","Line.Length2", "Center1", "Center2", "Start1","End1","Start2", "End2", "Positions1","Positions2","MutsCount1","MutsCount2","MutsCountTotal", "Z.Score", "Within.Range1", "Within.Range2","Intersection")
    }
    
    
    return(results.expanded)
  }
}
