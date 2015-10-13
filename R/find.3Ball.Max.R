find.3Ball.Max <-
function(sample, type = "first"){
  
  total.found = 0
  
  sample$MutsCount <- unlist(sample$MutsCount)
  sample <- sample[order(sample$MutsCount,decreasing=TRUE),]
  v<- sample$MutsCount
  l<- length(v)
  
  max.matrix <- matrix(data = 0, nrow = 0, ncol = 4)
  current.max <- 0
  colnames(max.matrix)<- c("i","j", "k", "sum")
  
  starti = 3
  startj = 2
  startk = 1
  
  cand<- as.data.frame(matrix(data = 0, nrow = 0, ncol =4))
  cand <- rbind(cand, c(v[starti] + v[startj] + v[startk],starti,startj,startk), deparse.level =0)
  names(cand)<- c("sum","i","j", "k")
  
  while(dim(cand)[1]!= 0){
    index <- which.max(cand$sum)
    sum  <- cand$sum[index]
    i <- cand$i[index]
    j <- cand$j[index]
    k <- cand$k[index]
    cand <- cand[-index,]
    
    intersectionIJ <- intersect(unlist(sample$Within.Range[i]), unlist(sample$Within.Range[j]))
    intersectionIK <- intersect(unlist(sample$Within.Range[i]), unlist(sample$Within.Range[k]))
    intersectionJK <- intersect(unlist(sample$Within.Range[j]), unlist(sample$Within.Range[k]))
    if ( ! ( (i==j) || (j==k) || (k==i)) & (length(intersectionIJ)==0) & (length(intersectionIK)==0) & (length(intersectionJK)==0)   ){
      #cat(sprintf("(%d,%d,%d,%d)\n", sum, i, j, k))
      total.found <- total.found +1
      if(total.found ==1){
        current.max <- sample$MutsCount[i]+sample$MutsCount[j] + sample$MutsCount[k]
      }
      
      if(type == "first"){
        break
      }else{
        if(sample$MutsCount[i]+sample$MutsCount[j] +sample$MutsCount[k]== current.max){ #still finding equivalent size balls
          max.matrix <- rbind(max.matrix, c(sample$Center[i],sample$Center[j], sample$Center[k],sample$MutsCount[i]+sample$MutsCount[j] +sample$MutsCount[k]))
          total.found <- total.found +1
        }
        else{
          break
        }
      }
    }
    
    if(k==1 && j==2  && i< l){
        cand <- rbind(cand, c( v[i+1]+v[j]+v[k],i+1, j, k ), deparse.level = 0)
        names(cand)<- c("sum","i","j", "k")
    }
    if(k==1 && j+1<i){
        cand <- rbind(cand, c( v[i]+v[j+1]+v[k],i, j+1, k), deparse.level = 0)
        names(cand)<- c("sum","i","j", "k")
    }
    if(k+1<j)
        cand <- rbind(cand, c(v[i]+v[j]+v[k+1], i, j ,k+1), deparse.level = 0)
        names(cand)<- c("sum","i","j", "k")
  }
  
  if(total.found == 0){ #This is an error case, all 2 balls overlap
    max.covered = -1
  }else{ #A success was found
    max.covered = sample$MutsCount[i]+sample$MutsCount[j]+sample$MutsCount[k]
  }
  
  if(type == "first"){
    return(max.covered)
  }else{
    
    if(max.covered != -1){
      results.expanded <- matrix(data = 0, nrow = dim(max.matrix)[1], ncol = 24)
      results.expanded <- as.data.frame(results.expanded)
      colnames(results.expanded) <- c("Line.Length1","Line.Length2", "Line.Length3", "Center1", "Center2","Center3", "Start1","End1","Start2", "End2","Start3","End3", "Positions1","Positions2","Positions3", "MutsCount1","MutsCount2", "MutsCount3","MutsCountTotal", "Z.Score", "Within.Range1", "Within.Range2", "Within.Range3", "Intersection")
      results.expanded$Center1 <- max.matrix[,1]
      results.expanded$Center2 <- max.matrix[,2]
      results.expanded$Center3 <- max.matrix[,3]
      results.expanded$MutsCountTotal <- max.matrix[,4]
        
      #fills in the matrix
      part1 <- lapply(results.expanded$Center1, function(x){sample[which(x == sample$Center) ,c("Line.Length","Start","End","Positions","MutsCount","Within.Range")]})
      part2 <- lapply(results.expanded$Center2, function(x){sample[which(x == sample$Center) ,c("Line.Length","Start","End","Positions","MutsCount","Within.Range")]})
      part3 <- lapply(results.expanded$Center3, function(x){sample[which(x == sample$Center) ,c("Line.Length","Start","End","Positions","MutsCount","Within.Range")]})
      results.expanded[,c("Line.Length1", "Start1","End1", "Positions1", "MutsCount1","Within.Range1")] <-  as.data.frame(t(sapply(part1,function(x){x})))
      results.expanded[,c("Line.Length2", "Start2","End2", "Positions2", "MutsCount2","Within.Range2")] <-  as.data.frame(t(sapply(part2,function(x){x})))
      results.expanded[,c("Line.Length3", "Start3","End3", "Positions3", "MutsCount3","Within.Range3")] <-  as.data.frame(t(sapply(part3,function(x){x})))
       
      #finds the ones with no intersection
      results.expanded$Intersection <- lapply(1:dim(results.expanded)[1], function(x){
      
        intersect12 <- intersect(unlist(results.expanded$Within.Range1[x]), unlist(results.expanded$Within.Range2[x]))
        intersect13 <- intersect(unlist(results.expanded$Within.Range1[x]), unlist(results.expanded$Within.Range3[x]))
        intersect23 <- intersect(unlist(results.expanded$Within.Range2[x]), unlist(results.expanded$Within.Range3[x]))
        return( union(union(intersect12,intersect13),intersect23))
      })
      no.intersect.indices <- which( unlist(lapply(results.expanded$Intersection, function(x){length(x)}))==0)
      results.expanded <- results.expanded[no.intersect.indices,]    
    }else{#if all 2 balls overlap, return a matrix with 0 rows which will be tested for later.
      results.expanded <- matrix(data = 0, nrow = 0, ncol = 24)
      results.expanded <- as.data.frame(results.expanded)
      colnames(results.expanded) <- c("Line.Length1","Line.Length2", "Line.Length3", "Center1", "Center2","Center3", "Start1","End1","Start2", "End2","Start3","End3", "Positions1","Positions2","Positions3", "MutsCount1","MutsCount2", "MutsCount3","MutsCountTotal", "Z.Score", "Within.Range1", "Within.Range2", "Within.Range3", "Intersection")
      
    }
    return(results.expanded)
  }
}
