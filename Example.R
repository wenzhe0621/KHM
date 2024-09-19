##====================================climate model======================================
# Source files 
source('Kernel PCA functions.R')
source('KPCA HM FUNCTIONS.R')
source("ErrorBar.R")
source("impLayoutplot.R")
source("HistoryMatching.R")
#=================================================================================================
# Load data 
#=================================================================================================
# Observation 
load("SANDU_LES_WAVE1.RData")
# Wave1 ensemble
load("SANDU_SCM_WAVE1.RData")
load("Wave1.RData")
# Get vertical grid
load("zgrid.RData")

data<-qvdata1<-qvdata
sample<-wave1sample<-wave_param_US[,2:6]
obs=lesdata[,,1]
obsref=lesdata[,,2]
timeN=72
nlevout=50
n=90

#Interpolate LES/SCM data on the chosen vertical grid (we use 1-72 times)
dim(data) <- c(timeN*nlevout, n)
dim(obs) <- c(timeN*nlevout, 1)
dim(obsref) <- c(timeN*nlevout, 1)

# Observation plots
clouldscol1<-rainbow(n=30,start=0.5001, end=1)
clouldscol2<-rainbow(n=15,start=0.1, end=0.4)
clouldscol<-c("darkgrey",clouldscol1,clouldscol2)

require('fields')
par(mfrow = c(1,1), mar=c(4,4,1,1))
image.plot(1:timeN,zgrid, matrix(t(obs),ncol=nlevout), ylim=c(0,max(zgrid )),col =clouldscol,ylab="Height (meters)",xlab="Time (hours)")

# Wave 1 emsemble plots
par(mfrow = c(5,9), mar=c(2,2.2,1,2))
for (i in 1:n) {
  image.plot(1:timeN,zgrid, matrix(t(data[,i]),ncol=nlevout), ylim=c(0,max(zgrid )),col =clouldscol,ylab="")
}

#========================================================================================================================
# expert's choice 
data_accept1=data_accept=c(4,6,7,12,17,19,23,30,37,43,44,71,77,83)
# plot the acceptable data
par(mfrow = c(3,5), mar=c(3,3,1,1))
for (i in data_accept) {
  image.plot(1:timeN,zgrid, matrix(t(data[,i]),ncol=nlevout), ylim=c(0,max(zgrid )),col =clouldscol,ylab="Height (meters)",xlab="Time (hours)")
}

# The cloud fraction for the later hour (time=68) of the simulation plot
obstime68<-lesdata[68,,1]
datatime68<-qvdata[68,,]
par(mfrow = c(1,1), mar=c(4,4,2,2))
plot(x=seq(1,length(obstime68),1),y=obstime68,type='l',lwd=2,col=2,xlab="Model level",ylab="Cloud fraction",ylim=c(min(data),max(c(data,obs))))
for(j in 1:n){
  lines(x=seq(1,length(obstime68),1),y=datatime68[,j],col=8,lwd=0.5)
}
for(j in data_accept){
  lines(x=seq(1,length(obstime68),1),y=datatime68[,j],col=3,lwd=2)
}
lines(x=seq(1,length(obstime68),1),y=obstime68,type='l',lwd=2,col=2)

#========================================================================================================================
# Plot for the first 14 ‘closest’ runs to the observation in model output space
obsE<-disE<-c()
for ( k in 1:length(obs)) {
  obsE[k] =sd(c(obs[k],obsref[k]))^2
}
dis<-sapply(1:n,function(i) t(obs-data[,i])%*%diag(1/(obsE+1))%*%(obs-data[,i]))
order(dis)[1:14]

par(mfrow = c(1,1), mar=c(2,4,2,2))
plot(x=seq(1,length(obstime68),1),y=obstime68,type='l',lwd=2,col=2,xlab="Model level",ylab="Cloud fraction",ylim=c(min(data),max(c(data,obs))))
for(j in 1:n){
  lines(x=seq(1,length(obstime68),1),y=datatime68[,j],col=8,lwd=0.5)
}
for(j in order(dis)[1:14]){
  lines(x=seq(1,length(obstime68),1),y=datatime68[,j],col=4,lwd=2)
}
lines(x=seq(1,length(obstime68),1),y=obstime68,type='l',lwd=2,col=2)

#========================================================================================================================
# kernel selection
opt1 <- GenSA(as.numeric(c(0.08011,0.4898,0.3943,0.2550)),optimdis,lower =c(0,0,0.00001,0.00001), upper = c(1,0.5,1,1),
             data_ini=data,obs_ini=obs, data_accept=data_accept,errvar=TRUE,control=list(threshold.stop=0.1)) 
#========================================================================================================================
# kernel pca
#par=opt1$par
par=c(0.08011,0.4898,0.3943,0.2550)
errvar <- calcA(expand.grid(1:72,1:50), d = c(par[3],par[4]), nu = 0)
Weight=diag(obsE)+ errvar 
Weight_inverse=chol2inv(chol(Weight))
Q_W=chol(Weight_inverse)

data_new<-Q_W%*%data
obs_new<-Q_W%*%obs

scale=max(Kernelmatrix_K(data_new,data_new,linearkernel))

k1<-function (x,y){
  (1-par[1])*scale*(exp(par[2]*(2*crossprod(x, y) -crossprod(x) -crossprod(y))))+par[1]*crossprod(x, y)}

kpca<-kernelpca(k1,data,n)
components<-kpca$C
V<-kpca$V
D<-kpca$D
(cumsum(D)/sum(D))
lamda1=(cumsum(D)/sum(D))>0.8
r1=min(which(lamda1 == TRUE))
r=min(r1,5)
obsprojection<-projectionXnew(as.vector(obs_new),data_new,k1,V)
acceptable_plots(sample,t(components)[,1:90],data_accept,obsprojection,r,main,xname=colnames(sample),yname=c("c1","c2","c3","c4","c5"))

KHM_pre1<-KHM_pre(par,data_ini=data,obs_ini=obs,data_accept=data_accept,errvar=TRUE)
acceptable_plots(sample,t(KHM_pre1$components)[,1:90],data_accept,KHM_pre1$obsprojection,KHM_pre1$r,main,xname=colnames(sample),yname=c("c1","c2","c3","c4","c5"))
KHM_pre1$threshold

dis<-sapply(1:n, function(k) discomponents(components[k,1:r],obsprojection[1:r]))
par(mfrow = c(1,1), mar=c(2,4,2,2))
hist(dis,breaks = 100)
optimf=c()
for (i in 1:n) {
  threshold=dis[i]
  optimf[i]=optimfunction(threshold, dis, data_accept)
}
max(optimf)
threshold=dis[order(-optimf)[1]]

# wave 1
# emulators
newN=1000
Newsample <- matrix(runif(n=5*newN,min=-1,max=1),ncol=5) 
colnames(Newsample) <- colnames(sample)
out<-list()
for (i in 1:r) {
  out[[i]] <-btgp(X=sample, Z=components[,i],XX=Newsample)
}
#leave one out plot
cname<-c('c1','c2','c3','c4','c5')

par(mfrow = c(r,5), mar = c(4,4,0.1,0.1))

par(mfrow = c(r,5),oma=c(0,0,0.1,0.5),mar = c(4,4,0,0),cex.lab=0.9, cex.axis=0.8)
for (k in 1:r) {
  for (i in 1:5) {
    PlotErrorBarPlot_inputsname(out[[k]]$Zp.mean,out[[k]]$Zp.s2,sample, components[,k] ,axis=i,observation=obsprojection[k],E1=0,E2=0,theta=2,ylab=cname[k],inputs=colnames(sample))
  }
}

# leave one out prediction
var0<-mean0<-list()
for (k in 1:n) {
  var0[[k]]<-diag(c(out[[1]]$Zp.s2[k],out[[2]]$Zp.s2[k],out[[3]]$Zp.s2[k],out[[4]]$Zp.s2[k],out[[5]]$Zp.s2[k]))
  mean0[[k]]<-c(out[[1]]$Zp.mean[k],out[[2]]$Zp.mean[k],out[[3]]$Zp.mean[k],out[[4]]$Zp.mean[k],out[[5]]$Zp.mean[k])
}
# Emulator prediction for newsampling
VAR<-MEAN<-c()
for (k in 1:newN) {
  VAR[[k]]<-diag(c(out[[1]]$ZZ.s2[k],out[[2]]$ZZ.s2[k],out[[3]]$ZZ.s2[k],out[[4]]$ZZ.s2[k],out[[5]]$ZZ.s2[k]))
  MEAN[[k]]<-c(out[[1]]$ZZ.mean[k],out[[2]]$ZZ.mean[k],out[[3]]$ZZ.mean[k],out[[4]]$ZZ.mean[k],out[[5]]$ZZ.mean[k])
}

impl<-c()
for (i in 1:newN){
  cutmean<-sum(VAR[[i]])
  cutvar<-2*sum((VAR[[i]])^2)+4*threshold*sum((VAR[[i]]))
  
  cut<-cutmean+3*sqrt(cutvar)+threshold
  proj.output <- MEAN[[i]]
  y <- as.vector(obsprojection[1:r] - proj.output)
  impl[i] <- crossprod(y,y)/cut
}
Cutoff<-1
par(mfrow = c(1,1), mar = c(4,4, 1, 1))
length(impl[impl<Cutoff])


impl<-c()
for (i in 1:newN){
  cutmean<-sum(VAR[[i]])
  cutvar<-2*sum((VAR[[i]])^2)
  
  cut<-cutmean+3*sqrt(cutvar)+threshold
  proj.output <- MEAN[[i]]
  y <- as.vector(obsprojection[1:r] - proj.output)
  impl[i] <- crossprod(y,y)/cut
}
Cutoff<-1
par(mfrow = c(1,1), mar = c(4,4, 1, 1))
length(impl[impl<Cutoff])


ImpData<-ImpData1<-cbind(Newsample,impl)
colnames(ImpData1)<-c(colnames(Newsample),"I")
VarNames <- colnames(Newsample)
breakVec <- sort(c(2000,seq(from=0,to=1,len=7)))
ImpList <- CreateImpList(whichVars = 1:6, VarNames = VarNames, ImpData = ImpData, Resolution=c(10,10), whichMax=1,Cutoff =Cutoff)
imp.layoutm11(ImpList,VarNames,VariableDensity=FALSE,newPDF=FALSE,the.title="",newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,tol=0.05)
length(ImpData1[,6][ImpData1[,6]<Cutoff])


#-------------------------------------
IMPL0<-c()
for (i in 1:n){
  V=sqrt(var0[[i]]+Id(r))
  y1 <- as.vector(obsprojection[1:r] - mean0[[i]])
  y<-backsolve(V,y1)
  IMPL0[i] <- crossprod(y,y)
}
max(IMPL0[data_accept])

IMPL2<-c()
for (i in 1:newN){
  V=sqrt(VAR[[i]]+Id(r))
  y1 <- as.vector(obsprojection[1:r] - MEAN[[i]])
  y<-backsolve(V,y1)
  IMPL2[i] <- crossprod(y,y)
}
Cutoff=max(IMPL0[data_accept])
length(IMPL2[IMPL2<Cutoff])

par(mfrow = c(1,1), mar = c(4,4, 1, 1))
ImpData<-ImpData1<-cbind(Newsample,IMPL2)
colnames(ImpData1)<-c(colnames(Newsample),"I")
VarNames <- colnames(Newsample)
breakVec <- sort(c(2000,seq(from=0,to=Cutoff,len=7)))
ImpList <- CreateImpList(whichVars = 1:6, VarNames = VarNames, ImpData = ImpData1, Resolution=c(10,10), whichMax=1,Cutoff =Cutoff)
imp.layoutm11(ImpList,VarNames,VariableDensity=FALSE,newPDF=FALSE,the.title="",newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,tol=0.5)

#save(ImpData1,file="ImpData1.RData")
