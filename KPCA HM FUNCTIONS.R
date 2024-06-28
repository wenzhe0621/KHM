packages <- c('tgp','kernlab','GenSA', 'parallel', 'cvTools', 'fields', 'mco', 'far', 'lhs', 'tensor', 'matrixcalc', 'rstan','MASS','scatterplot3d')
sapply(packages, require, character.only = TRUE)
#--------------------mutipred
pred.em <- function(em, new, newy,...){
  require(fields)
  x = em$x
  n = em$n
  regmodel = em$regmodel
  nu = em$nu
  error = em$error
  H = em$H
  q = em$q
  sigma2 = em$sigma2
  type <- em$type
  active <- em$active
  if (length(active) == 1){
    x.ac <- as.data.frame(x[,active])
    names(x.ac) <- active
  }
  else {
    x.ac <- x[,active]
  }
  if (is.null(dim(new)) == T){
    new <- new
    new.ac <- as.data.frame(new)
    colnames(new.ac) <- active
  }
  else if ((dim(new)[1]*dim(new)[2]) == 1){
    new <- new
    new.ac <- as.data.frame(new)
    colnames(new.ac) <- active
  }
  else if ((dim(new)[1]*dim(new)[2]) == dim(new)[2] & length(active) == 1){
    new.ac <- as.data.frame(new[,active])
    names(new.ac) <- active
  }
  else {
    new = new[,names(x)]
    new.ac <- as.data.frame(new[,active])
    if (length(active) == 1){
      names(new.ac) <- active
    }
  }
  if (type == "BR"){
    if (nu == 1){
      postm = predict(regmodel,new.ac)
      tt = terms(regmodel)
      Terms = delete.response(tt)
      mm = model.frame(Terms, new.ac, xlev = regmodel$xlevels)
      hx = model.matrix(Terms, mm, contrasts.arg = regmodel$contrasts)
      postcov = sigma2 * (diag(dim(new.ac)[1]) + (hx)%*%chol2inv(
        chol(t(H)%*%H))%*%t(hx))
    }
    else {
      d = em$d
      tx = calctx(x.ac,new.ac,d,nu,...)
      Q = em$Q
      x.c = em$x.c
      h.c = em$h.c
      comp2 = em$comp2
      y.c = backsolve(Q, tx, transpose = TRUE)
      mean.adj = crossprod(y.c,x.c)
      cov.adj = crossprod(y.c,y.c)
      priormean = pred.gls(regmodel, new.ac)
      priorcov = calcA(new.ac, d, nu,...)
      postm = priormean + mean.adj
      postcov = priorcov - cov.adj
      tt = terms(regmodel)
      Terms = delete.response(tt)
      mm = model.frame(Terms, new.ac, xlev = regmodel$xlevels)
      hx = model.matrix(Terms, mm, contrasts.arg = regmodel$contrasts)
      comp1 = hx - crossprod(y.c, h.c)
      postcov = sigma2 * (postcov + comp1 %*% comp2 %*% t(comp1))
    } }
  return(list(post.m = postm, post.cov = diag(postcov)))
}

MultiPred <- function(ems, design){
  Expectation <- matrix(0, nrow = dim(design)[1], ncol = length(ems))
  Variance <- matrix(0, nrow = dim(design)[1], ncol = length(ems))
  s <- dim(design)[1]/10
  for (i in 1:s){
    EmOutput <- lapply(1:length(ems), function(e) pred.em(ems[[e]],design[(10*(i-1) + 1):(10*i),]))
    for (j in 1:length(ems)){
      Expectation[(10*(i-1) + 1):(10*i),j] <- EmOutput[[j]]$post.m
      Variance[(10*(i-1) + 1):(10*i),j] <- EmOutput[[j]]$post.cov
    } }
  return(list(Expectation=Expectation, Variance=Variance))
}


#-------------------------Next wave sampling --------------------------------------------------------------------------------
#This function only works for moving signal example.
NewWAVE_SAMPLING<-function(X.NEW,N,data,k1,V){
  designwaveNEW<-sample(nrow(X.NEW),N,rep=F)
  waveNEWsample<-as.data.frame(X.NEW[designwaveNEW,])
  
  waveNEWData<-matrix(0,ncol = N,nrow = 100)
  waveNEWcol<- vector()
  for (i in 1:n){
    waveNEWData[,i] <- fn(as.numeric((waveNEWsample[i,]+1)/2))$fx
    waveNEWcol[i] <- fn(as.numeric((waveNEWsample[i,]+1)/2))$shape
  }
  waveNEWcomponents<-t(sapply(1:50,function(i) projectionXnew(waveNEWData[,i],data,k1,V)))
  
  return(list(waveNEWsample=waveNEWsample,waveNEWData=waveNEWData,waveNEWcol=waveNEWcol,waveNEWcomponents=waveNEWcomponents))
}

#-------------------------Emulator tgp--------------------------------------------------------------------------------
tgpemulator<-function(wavenewsample,wavenewcomponents,X.new,obsprojection,r){
  wavenewout<-list()
  for (i in 1:r) {
    wavenewout[[i]] <-btgp(X=wavenewsample, Z=wavenewcomponents[,i],XX=X.new)
  }
  return(wavenewout)
}

#-------------------------DISCREPANCY AND ERROR--------------------------------------------------------------------------------
Error_var<-function(errvar,n,r){
  errvar <- calcA(expand.grid(1:10,1:10), d = c(1,1), nu = 0)
  esample <- matrix(mvrnorm(n, c(rep(0,n)), errvar),nrow=n)
 
  Keprojection<-sapply(1:n ,function(i) projectionXnew(esample[,i],data,k1,V))
  colnames(Keprojection)<-paste("e",1:n,sep="")
  VARe<-diag(cov(t(Keprojection[1:r,])))
  W<-diag(VARe)
  
  return(var<-W)
}

#------------------------Calculating implausibilitie -------------------------------------------------------------
###  mahalanobis distance
implCoe<-function(Var,mean,obsprojection,newN,Error,Disc){
  if(missing(Error) & missing(Disc)){
  impl <-c()
  for (i in 1:newN){
    Var <- var[[i]]
    Q <- sqrt(Var)
    proj.output <- mean[[i]]
    y <- backsolve(Q, as.vector(obsprojection[1:r] - proj.output), transpose = TRUE)
    impl[i] <- crossprod(y,y)
  }
  }else{
    impl <-c()
    for (i in 1:newN){
      Var <- var[[i]]+Error+Disc
      Q <- sqrt(Var)
      proj.output <- mean[[i]]
      y <- backsolve(Q, as.vector(obsprojection[1:r] - proj.output), transpose = TRUE)
      impl[i] <- crossprod(y,y)
    }
  }
  return(impl)
}


###  L2 norm distance
implL2Coe<-function(obsprojection,mean,var,r,newN,b){
  dist1<-E<-Var<-c()
  for (i in 1:newN){
    dist1[i]<-t(obsprojection[1:r]-mean[[i]])%*%(obsprojection[1:r]-mean[[i]])
    E[i]<-b+(1+2*b)*sum(diag(var[[i]]))
    Var[i]<-(1+2*b)^2*sum(diag(var[[i]]))
  }
  threshold<-E+3*Var
  impl<- dist1/threshold
  return(list(dist=dist1,impl=impl))
}

#-----------------------------------------------------------------------------------------------------------------
HM_hmwave1<-function(impl,Newsample,Cutoff){
  ImpData<-cbind(Newsample,impl)
  colnames(ImpData)<-c(names(Newsample),"I")
  VarNames <- names(Newsample)
  ImpList <- CreateImpList(whichVars = 1:5, VarNames = VarNames, ImpData = ImpData, Resolution=c(10,10), whichMax=1,Cutoff = Cutoff)
  imp.layoutm11(ImpList,VarNames,VariableDensity=FALSE,newPDF=FALSE,the.title="",newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE)
  
  NROY1<- impl<=Cutoff
  X.2<-Newsample[NROY1,]
  newN2<-dim(X.2)[1]
  return(list(ImpData=ImpData,X.next=X.2,NROYlength=newN2))
}


NEW_hmwaveM<-function(implM,Newsample,Cutoff,ImpDatawave1){
  # n<-dim(Newsample)[2]
  
  ImpData<-ImpDatawave1
  ImpData[,6][which(ImpData[,6]<=Cutoff)]<-implM
  
  #length(ImpData[,n+1][which(ImpData[,n+1]<=Cutoff)])
  #length(wave2$impl)
  
  colnames(ImpData)<-c(names(Newsample),"I")
  VarNames <- names(Newsample)
  
  ImpList <- CreateImpList(whichVars = 1:5, VarNames = VarNames, ImpData = ImpData, Resolution=c(10,10), whichMax=1,Cutoff = Cutoff)
  imp.layoutm11(ImpList,VarNames,VariableDensity=FALSE,newPDF=FALSE,the.title="",newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE)
  
  NROY<- implM< Cutoff
  X.new<-Newsample[NROY,]
  newN<-dim(X.new)[1]
  return(list(ImpData=ImpData,X.next=X.new,NROYlength=newN))
}



calcA <- function(x,d,nu,cov = "gauss",d2 = c(rep(2,length(d)))){
  n = dim(x)[1] # number of training points
  p = dim(x)[2] # number of parameters
  if (is.null(n) == TRUE){
    n <- length(x)
  }
  if (is.null(p) == TRUE){
    p <- 1
  }
  stopifnot(p == length(d))
  if (p == 1){
    x <- x }
  else {
    for (i in 1:p){
      x[,i] = as.numeric(as.character(x[,i]))
    }
  }
  if (cov == "gauss") {
    leng = as.matrix(dist(scale(x,center=FALSE,scale=d),method="euclidean",
                          diag=TRUE,upper=TRUE))
    A = nu*diag(n) + (1 - nu)*exp(-(leng^2))
  }
  return(A) }
