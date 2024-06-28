PlotErrorBarPlot<- function(predict,V, x, y ,axis=3, heading='',observation,a,b,theta,E1,E2,ylab='') {
  # Plots the emulator predictions together with error bars and true values.
  #
  # Args: 
  #  prediction: by the emulation.
  #  X: A normalized validation matrix.
  #  y: A toy function evaluated at X.
  #
  # Returns:
  #  A plot of emulator predicitons together with error bars and true values
  if (missing(a)) {
    a=-1
  }
  if (missing(b)) {
    b=1
  }
  if (missing(E1)) {
    E1=0.01
  }
  if (missing(E2)) {
    E2=0.05
  }
  if (missing( ylab)) {
    ylab='f(x)'
  }

  
  N <- dim(x)[1]
  outside <- c()
  inside <- c()
  stand.dev<-sqrt(V)
  maxval <- predict+theta*stand.dev
  minval <- predict-theta*stand.dev
  minv <-min (minval)
  maxv<-max(maxval,observation+theta*(E1+E2)^(1/2))
  
  for(i in 1:length(y)) {
    if((minval[i]) <= y[i] & y[i] <= (maxval[i])){
      inside <- c(inside,i)
    }
  }
  outside <- c(1:N)[-inside]
  
  if (missing(observation)){
    plot(x[ , axis], predict, pch=20, ylab=ylab, xlab=paste('x', toString(axis), sep=""), ylim=c(minv, maxv),
         xlim=c(a,b), main=heading)
    
    for(i in 1:N){
      arrows(x[i, axis],minval[i] , x[i, axis],maxval[i],
             length=0.05, angle=90, code=3, col='black')
    }
    points(x[outside, axis], y[outside], pch=20, col='red')
    points(x[inside, axis], y[inside], pch=20, col='deepskyblue3')
    
  } else {
    
    plot(x[ , axis], predict, pch=20, ylab=ylab, xlab=paste('x', toString(axis), sep=""), ylim=c(minv, maxv),
         xlim=c(a,b), main=heading)
    
    for(i in 1:N){
      arrows(x[i, axis],minval[i] , x[i, axis],maxval[i],
             length=0.05, angle=90, code=3, col='black')
    }
    points(x[outside, axis], y[outside], pch=20, col='red')
    points(x[inside, axis], y[inside], pch=20, col='deepskyblue3')
    
    abline(h= as.numeric(observation),col='pink')
    abline(h= as.numeric(observation) + theta*(E2+E1)^(1/2),col='pink',lty=2)
    abline(h= as.numeric(observation) - theta*(E1+E2)^(1/2),col='pink',lty=2)
    
  }
}









