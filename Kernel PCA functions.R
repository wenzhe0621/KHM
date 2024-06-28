#kernel package
require("kernlab")

### manualkernel functions ###

#  gaussxian kernel
kernelGaussian =function (x,y) exp(-crossprod(data[,1]-obs)/1000)

# Polynomial Kernel
KernelPolynomial = function (x,y)  (t(x) %*% y)^2

# Sigmoid Kernel
KernelSigmoid<- function(x,y) (t(x) %*% y + 1)

# Inverse Multiquadric Kernel
KernelInverse<- function(x,y) 1/sqrt(norm(matrix(x-y), type="F") + 1)

# degree-2 polynomials kernel
KernelPolynomial_2 <- function (x,y)  (t(x) %*% y + 1)^2
# degree-2 polynomials feature map function
pi<- function(x){
  x<-matrix(x,ncol=1)
  k<-nrow(x)
  c1<- sapply(1:k, function(k) x[k,]^2)
  
  c2<-matrix(0,k,k) 
  for (i in 1:k) {
    for (j in 1:k) {
      c2[i,j]<-sqrt(2)*x[i,]*x[j,]
    }
  }
  
  c3<- sapply(1:k, function(k) sqrt(2)*x[k,])
  fi<-c(c1,c2[lower.tri(c2,diag = FALSE)],c3,1)
  return(fi)
}

###  Kernel PCA  ###
kernelpca<-function(k1,data,r){
  # k1: kernel function
  # data: simulation outputs data
  # r: the first r coefficients 
  # RETURN: the projection，eigenvector and eigenvlues
  n<-ncol(data)
  if (missing(r)){r=n}
  K<-Kernelmatrix_K(data,data,k1)
  
  ones <-1/n* matrix(1, n, n)
  H <- Id(n)- ones
  K_norm<- H%*%K%*%H 
  #K_norm <- K-ones%*%K-K%*%ones+ones%*%K%*%ones
  
  #caculated the eigenvectors and eigenvalues of centred kernel matrix
  res<-eigen(K_norm)
  V<-res[["vectors"]]
  D<-res[["values"]]
  
  #caculated covariance erigenvectore and eigenvalues
  D_norm<-D
  V_norm<-V
  for (i in 1:r) {
    if (D[i] < 1e-6) { break }
    V_norm[,i] = V[,i] / sqrt(D_norm[i])
  }
  C<- K_norm %*% V_norm[,1:r]
  return(list(C=C,V=V_norm[,1:r],D=D_norm[1:r],K=K,K_norm=K_norm))
}

### Kernel matrix ###
Kernelmatrix_K<-function(data1,data2,k1){
  #note: dim(data)= (data dimension, data numeber)
  #n:  dim of data1 
  #m:  dim of data2
  #k1:  kernel function 
  
  n<-ncol(data1)
  m<-ncol(data2)
  
  if (is.null(ncol(data1)) == T & is.null(ncol(data2)) == F){
    K<-matrix(0, ncol = 1, nrow = m)
    for (i in 1:m) {
      K[i,] = k1(data1, data2[,i])
    }
  }
  else if(is.null(ncol(data2)) == T &is.null(ncol(data1)) == F){
    K<-matrix(0, ncol = 1, nrow = n)
    for (i in 1:n) {
      K[i,] = k1(data1[,i], data2)
    }
  }
  else if(is.null(ncol(data2)) == T &is.null(ncol(data1)) == T){
    K<- k1(data1, data2)
  }
  else{
    K = matrix(0, ncol = n, nrow = m)
    for (i in 1:n) {
      for (j in 1:m) {
        K[i,j] = k1(data1[,i], data2[,j])
      }}
  }
  return(K)
}

### centred Kernel function ###
kernel_norm<-function(data,k1,x1,x2){
  
  if (missing(x2)) {
    n<-ncol(data)
    one <-  matrix(1, n, 1)
    ones <- 1/n* matrix(1, n, n)
    H <- Id(n)- ones
    
    Kx1<-Kernelmatrix_K(as.vector(x1),data,k1)
    K<-Kernelmatrix_K(data,data,k1)
    
    k<-H%*%(Kx1-1/n*K%*%one)
    
  }else{
    n<-ncol(data)
    one<-  matrix(1, n, 1)
    Kx1<-Kernelmatrix_K(as.vector(x1),data,k1)
    Kx2<-Kernelmatrix_K(as.vector(x2),data,k1)
    K<-Kernelmatrix_K(data,data,k1)
    
    k<-k1(x1,x2)-
      1/n*t(one)%*%Kx1-
      1/n*t(one)%*%Kx2+
      1/n^2*t(one)%*%K%*%one
  }
  return(k)
}

### Indentity function ###
Id<-function(n){
  diag(c(1),nrow=n,ncol=n)
}

### feature space distance with known f(x) ### 
featureDISTANCE<-function(x1,x2,k1){
  # x1,x2 are the output on feature space 
  # k1 is choosen  kernel
  term1<-k1(x1,x1)
  term2<-k1(x2,x2)
  term3<-2*k1(x1,x2)
  dis<-term1+term2-term3
  return(dis)
}

### feature space distance based on coefficients ###
featureDISTANCE_coe<-function(data,obs,C,k1,r){
  V_norm<-kernelpca(k1,data,r)$V
  
  term1<-kernel_norm(data,k1,obs,obs)
  term2<-t(C)%*%C
  term3<-C%*%t(V_norm)%*%kernel_norm(data,k1,obs)
  
  dis<-term1+term2-2*term3
  #return(dis)
  return(dis)
}

### plot field ###
plot.field <- function(field, dim1, col = rainbow(100,start=0.1,end=0.8),zlim, ...){
  require(fields)
  x <- 1:dim1
  y <- 1:dim1
  if (missing(zlim)){
    image.plot(x, y, matrix(field, nrow = dim1), col = col, add = FALSE, ...)
  }
  else{
    image.plot(x, y, matrix(field, nrow = dim1), col = col,zlim=zlim, add = FALSE, ...)
  }
}

### projection of observation function ###
projectionXnew<-function(x1,data,k1,V_norm){
  projection<-t(V_norm)%*%kernel_norm(data,k1,x1)
  return(projection)
}

### projection test function ###
Centredk_x1x2<-function(x1,x2,k1,data){
  n<-ncol(data)
  one<-  matrix(1, n, 1)
  k_x1x2<-Kernelmatrix_K(x1,x2,k1)-
    1/n*t(one)%*%Kernelmatrix_K(x2,data,k1)-
    1/n*t(one)%*%Kernelmatrix_K(x1,data,k1)+
    1/n^2*(t(one)%*%K%*%one)
  return(k_x1x2)
}

### PRE IMAGE FUNCTION ###
NEW_z<-function(pre_z,data,gamma,k1){
  term1<-matrix(0, ncol = 1, nrow = n)
  for (i in 1:n) {
    term1[i,] = gamma[i,]*k1(data[,i],pre_z)
  }
  z=(data%*% term1)/sum(term1)
  return(z)
}

pre_image<-function(data,rank,eigVector,y,k1){
  
  # Reconstruction fix-point interation function
  #   y: dimensionanlity-reduced data
  #	  eigVector: eigen-vector obtained in kPCA
  #   data: data matrix
  #   para: parameter of Gaussian kernel
  #  	z: pre-image of y
  ####TEST
  #eigVector<-V_norm
  
  N= ncol(data)
  datamean <- apply(data, 1, mean)
  pre_z<-datamean
  
  gamma=matrix(0,ncol =1,nrow = N)
  for (i in 1:n) {
    gamma[i,]=eigVector[i,1:rank]%*%y
  }
  
  z<-NEW_z(pre_z,data,gamma,k1)
  while(norm(pre_z-z)/norm(z)>0.00001){
    pre_z<-z
    z<-NEW_z(pre_z,data,gamma,k1)
    if (norm(pre_z-z)/norm(z)<0.00001){break}
  }
  return(z)
}

