
UniImplausibilityStan <- function(NewData, Emulator, Discrepancy, Obs, ObsErr, is.GP=NULL, FastVersion=FALSE, batches=500,multicore=FALSE){
  if(is.null(is.GP)){
    is.GP <- is.list(Emulator$ParameterSamples)
  }
  if(multicore){
    if(!is.GP){
      tEmulator <- EMULATOR.gpstan.multicore(NewData, Emulator, GP=FALSE, FastVersion = FastVersion,batches = batches)
    }
    else{
      tEmulator <- EMULATOR.gpstan.multicore(NewData, Emulator=Emulator, GP=TRUE, FastVersion = FastVersion,batches=batches)
    }
  }
  else{
    if(!is.GP){
      tEmulator <- EMULATOR.gpstan(NewData, Emulator, GP=FALSE, FastVersion = FastVersion)
    }
    else{
      tEmulator <- EMULATOR.gpstan(NewData, Emulator=Emulator, GP=TRUE, FastVersion = FastVersion)
    }
  }
  abs(Obs - tEmulator$Expectation)/sqrt(tEmulator$Variance + Discrepancy + ObsErr)
}


ManyImplausibilitiesStan <- function(NewData, EmulatorList, Discrepancy, Obs, ObsErr, is.GP=NULL, FastVersion=FALSE,multicore=FALSE,batches=500){
  if(!is.null(is.GP)){
    stopifnot(length(is.GP)==length(EmulatorList))
  }
  else{
    is.GP  <- unlist(lapply(EmulatorList, function(e) is.list(e$ParameterSamples)))
  }
  sapply(1:length(EmulatorList), function(k) UniImplausibilityStan(NewData, EmulatorList[[k]], Discrepancy[k], Obs[k], ObsErr[k], is.GP[k], FastVersion=FastVersion,multicore=multicore,batches=batches))
}

RuledOut <- function(ImpVec,whichMax=1, Cutoff=3){
  length(which(ImpVec>Cutoff)) >= whichMax
}

NROY <- function(Imp.Matrix, whichMax=3, Cutoff=3){
  sapply(1:length(Imp.Matrix[,1]), function(j) !RuledOut(Imp.Matrix[j,],whichMax,Cutoff))
}

MaxImp <- function(implausibilities,whichMax){
  if(whichMax==1)
    return(max(implausibilities))
  else{
    return(sort(implausibilities,decreasing=TRUE)[whichMax])
  }
}

ImpDensityPanel <- function(x1var, x2var, ImpData, nEms, Resolution=c(10,10), whichMax=3,Cutoff=3){
  #Find all points relevant to a pixel
  x1RHboundaries <- seq(from=-1,to=1,length=Resolution[1]+1)[-1]
  x2RHboundaries <- seq(from=-1,to=1,length=Resolution[2]+1)[-1]
  x1LHboundaries <- c(-1,x1RHboundaries[-Resolution[1]])
  x2LHboundaries <- c(-1,x2RHboundaries[-Resolution[2]])
  whereX1 <- which(colnames(ImpData)==x1var)
  whereX2 <- which(colnames(ImpData)==x2var)
  n <- prod(Resolution)
  xs <- 1:Resolution[1]
  ys <- 1:Resolution[2]
  xyGrid <- expand.grid(xs,ys)
  ImpOnePixel <- function(index){
    tx1Indices <- which(ImpData[,whereX1] > x1LHboundaries[xyGrid[index,1]] & ImpData[,whereX1] < x1RHboundaries[xyGrid[index,1]])
    tx2Indices <- which(ImpData[,whereX2] > x2LHboundaries[xyGrid[index,2]] & ImpData[,whereX2] < x2RHboundaries[xyGrid[index,2]])
    tIndices <- tx1Indices[which(tx1Indices %in% tx2Indices)]
    N <- length(tIndices)
    if(N<1){
      #No points in our design in this column
      return(c(NA,NA))
    }
    else{
      Timps <- ImpData[tIndices,seq(from=ncol(ImpData),by=-1,length.out = nEms)]
      if(!is.matrix(Timps))
        Timps <- as.matrix(Timps,ncol=1)
      tNROY <- NROY(Timps,whichMax,Cutoff)
      if(sum(tNROY)==0)
        tMaxs <- apply(Timps,1,MaxImp,whichMax=whichMax)
      else{
        TimpsFinal <- Timps[tNROY,]
        if(is.null(dim(TimpsFinal)))
          TimpsFinal <- matrix(TimpsFinal, nrow=1)
        tMaxs <- apply(TimpsFinal,1,MaxImp,whichMax=whichMax)
      }
      c(min(tMaxs),sum(tNROY)/N)
    }
  }
  VectOnePixel <- Vectorize(ImpOnePixel)
  VectOnePixel(1:n)
}

#Inputs:
#Unlike for wave 1, ImpData now has a list for each wave.
#List elements: 
#1. Design (this is the big N-point design in all input space for which we have implausibilities)
#2. NROY. A matrix N x M giving the NROY status at each wave. 
#E.g. Point k ruled out at wave 2 (not wave 1) has ImpData$NROY[k,] = c(TRUE, FALSE).
#Note c(FALSE, TRUE not possible)
#3. Impl. A list of M implausibility matrices
ImpDensityPanelWaveM <- function(x1var, x2var, ImpData, nEms, Resolution=c(10,10), whichMax=3){
  if(!is.list(ImpData)){
    stop("ImpData must be a list for multi-wave history matching, see code.")
  }
  #Find all points relevant to a pixel
  x1RHboundaries <- seq(from=-1,to=1,length=Resolution[1]+1)[-1]
  x2RHboundaries <- seq(from=-1,to=1,length=Resolution[2]+1)[-1]
  x1LHboundaries <- c(-1,x1RHboundaries[-Resolution[1]])
  x2LHboundaries <- c(-1,x2RHboundaries[-Resolution[2]])
  whereX1 <- which(colnames(ImpData$Design)==x1var)
  whereX2 <- which(colnames(ImpData$Design)==x2var)
  n <- prod(Resolution)
  xs <- 1:Resolution[1]
  ys <- 1:Resolution[2]
  xyGrid <- expand.grid(xs,ys)
  ImpOnePixel <- function(index){
    tx1Indices <- which(ImpData$Design[,whereX1] > x1LHboundaries[xyGrid[index,1]] & ImpData$Design[,whereX1] < x1RHboundaries[xyGrid[index,1]])
    tx2Indices <- which(ImpData$Design[,whereX2] > x2LHboundaries[xyGrid[index,2]] & ImpData$Design[,whereX2] < x2RHboundaries[xyGrid[index,2]])
    tIndices <- tx1Indices[which(tx1Indices %in% tx2Indices)]
    N <- length(tIndices)
    if(N<1){
      #No points in our design in this column
      return(c(NA,NA))
    }
    else{
      M <- length(ImpData$NROY[1,])
      tNROY <- ImpData$NROY[tIndices,M]
      Timps <- ImpData$Impl[[1]][tIndices,]
      if(!is.matrix(Timps))
        Timps <- as.matrix(Timps,ncol=1)
      #for(j in 2:M){ # serious bug in this bit of the code!!!
      #  nroyindices <- which(!ImpData$NROY[tIndices,j-1])
      #  Timps[nroyindices,] <- ImpData$Impl[[j-1]][tIndices[nroyindices],]
      #}
      for(j in 2:M) {
        nroyindices <- which(ImpData$NROY[tIndices,j-1])
        Timps[nroyindices, ] <- ImpData$Impl[[j]][tIndices[nroyindices], ]
      }
      if(sum(tNROY)==0)
        tMaxs <- apply(Timps,1,MaxImp,whichMax=whichMax)
      else{
        TimpsFinal <- Timps[tNROY,]
        if(is.null(dim(TimpsFinal)))
          TimpsFinal <- matrix(TimpsFinal, nrow=1)
        tMaxs <- apply(TimpsFinal,1,MaxImp,whichMax=whichMax)
      }
      c(min(tMaxs),sum(tNROY)/N)
    }
  }
  VectOnePixel <- Vectorize(ImpOnePixel)
  VectOnePixel(1:n)
}

#whichVars indexes all of the variables you want to plot implausibility for, IN ORDER
CreateImpList <- function(whichVars, VarNames, ImpData, nEms=1, Resolution=c(15,15), whichMax=3,Cutoff=3){
  combGrid <- expand.grid(whichVars[-length(whichVars)],whichVars[-1])
  badRows <- c()
  for(i in 1:length(combGrid[,1])){
    if(combGrid[i,1] >= combGrid[i,2])
      badRows <- c(badRows,i)
  }
  combGrid <- combGrid[-badRows,]
  combGrid <- combGrid[do.call(order,combGrid),]
  gridList <- lapply(whichVars[-length(whichVars)], function(k) combGrid[which(combGrid[,1]==k),])
  lapply(gridList, function(e) lapply(1:length(e[,1]), function(k) 
    ImpDensityPanel(x1var=VarNames[e[k,1]], x2var=VarNames[e[k,2]], 
                    ImpData = ImpData, nEms=nEms, Resolution=Resolution, 
                    whichMax = whichMax,Cutoff=Cutoff)))
}
# CreateImpList for wave M.
CreateImpListWaveM <- function(whichVars, VarNames, ImpData, nEms=1, Resolution=c(15,15), whichMax=3) {
  combGrid <- expand.grid(whichVars[-length(whichVars)],whichVars[-1])
  badRows <- c()
  for(i in 1:length(combGrid[,1])){
    if(combGrid[i,1] >= combGrid[i,2])
      badRows <- c(badRows,i)
  }
  combGrid <- combGrid[-badRows,]
  combGrid <- combGrid[do.call(order,combGrid),]
  gridList <- lapply(whichVars[-length(whichVars)], function(k) combGrid[which(combGrid[,1]==k),])
  lapply(gridList, function(e) lapply(1:length(e[,1]), function(k) 
    ImpDensityPanelWaveM(x1var=VarNames[e[k,1]], x2var=VarNames[e[k,2]], ImpData = ImpData, nEms=nEms, Resolution=Resolution, whichMax = whichMax)))
}
ImplField <- function(Basis, Expectation, Variance, Obs, Error, Disc){
  proj.output <- Expectation
  recon.output <- Recon(proj.output, Basis)
  var.output <- diag(Variance)
  recon.var <- Basis %*% var.output %*% t(Basis)
  V <- Error + Disc + recon.var
  Q <- chol(V)
  y <- backsolve(Q, as.vector(Obs - recon.output), transpose = TRUE)
  impl <- crossprod(y,y)
  return(impl)
}


BasisImplausibility <- function(Basis, EmPreds, Obs, Error, Disc){
  N <- length(EmPreds[[1]]$Expectation)
  getImpl <- function(index){
    tEx <- unlist(lapply(EmPreds, function(e) e$Expectation[index]))
    tVar <- unlist(lapply(EmPreds, function(e) e$Variance[index]))
    ImplField(Basis, tEx, tVar, Obs, Error, Disc)
  }
  Imps <- rep(0,N)
  gimp <- Vectorize(getImpl)
  Imps[1:N] <- gimp(1:N)
  Imps
}



# Function to map the NROY indices from each wave to the original space
MapToOroginalSpace <- function(NROY.list) {
  M <- length(NROY.list)
  if(M == 2) NROY.map.to.origin <- NROY.list[[M-1]][NROY.list[[M]]]
  else {
    NROY.map.to.origin <- NROY.list[[M-1]][NROY.list[[M]]]
    for(i in 2:(M-1)) NROY.map.to.origin = NROY.list[[M - i]][NROY.map.to.origin]
  }
  return(NROY.map.to.origin)
}
# Function to construct ImpData for ImpDensityPanelWaveM function
ImpDataWaveM <- function(Design, NROY.list, Impl.list) {
  origin.space <- 1:dim(Design)[1]
  M <- length(NROY.list) # number of waves
  NROYM <- matrix(rep(Impl.list[[1]]<3, M), ncol = M)
  Timps <- list()
  NROYL <- list()
  Timps[[1]] <- matrix(Impl.list[[1]], ncol = 1)
  NROYL[[1]] <- origin.space
  for(i in 2:M) {
    NROYL[[i]] <- NROY.list[[i-1]]
    NROYL.origin <- MapToOroginalSpace(NROYL)
    NROYM[NROYL.origin, i:M] <- Impl.list[[i]] < 3
    Timps[[i]] <- matrix(rep(NA, dim(Design)[1]), ncol = 1)
    Timps[[i]][NROYL.origin, 1] <- Impl.list[[i]]
  }
  return(list(Design = Design, NROY = NROYM, Impl = Timps))
}

