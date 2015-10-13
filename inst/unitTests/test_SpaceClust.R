test_SpaceClust  <- function() {
  Residue <- c("MET","THR","GLU","TYR")
  Can.Count<- c(1,2,3,4)
  SideChain <- c("A","A","A","A")
  XCoord <- c(62.935,63.155,65.289,64.899)
  YCoord <- c(97.579,95.525,96.895,96.220)
  ZCoord <- c(30.223,27.079,24.308,20.615)
  
  mutation.matrix <-as.data.frame(matrix(data = c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),ncol = 4))
  names(mutation.matrix)<- c("V1","V2","V3","V4")
  PositionList <- data.frame(Residue = Residue, Can.Count = Can.Count, SideChain = SideChain, XCoord = XCoord, YCoord = YCoord, ZCoord = ZCoord)
  
  resultBallSim <- SpaceClust(mutation.matrix, PositionList, numsims =1000, radii.vector = c(1,2,3,4) , alpha = .05, method = "SimMax")
  checkEquals(length(resultBallSim), 14)
  checkTrue(resultBallSim$optimal.num.spheres>0)
  checkTrue(resultBallSim$optimal.radius > 0)
  
  resultBallPoisson <- SpaceClust(mutation.matrix, PositionList, numsims =1000, radii.vector = c(1,2,3,4) , alpha = .1, method = "Poisson")
  checkEquals(length(resultBallPoisson), 2)
  checkTrue(resultBallPoisson$best.radii> 0)
}