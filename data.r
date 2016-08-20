rmultnorm<-function(n,muvec,sigmat){
  # n is the number of random vectors to be generated
  # mu is the mean vector, sigmat is the variance-covariance matrix
  # the function returns an n by p matrix with each row being a
  # random vector from a multivariate normal with mean vector muvec
  # and covariance matrix sigmat
  if(length(muvec)==1){
    temp<-rnorm(n,muvec,sqrt(sigmat))
    return(temp)
  }
  else{
    sigeigen<-eigen(sigmat)
    amat<-sigeigen$vectors%*%diag(sqrt(sigeigen$values))
    temp<-matrix(rnorm(n*length(muvec),0,1),ncol=n)
    temp<-(amat%*%temp)+muvec
    temp<-t(temp)
    return(temp)
  }
}


generate <- function(k)
{
  beta <- c(rep(0,96),c(0.1,0.1,0.1,0.1),c(1),c(0.6,-0.8,0.4),c(0.1,0.1,0.1,0.1),rep(0,3892),rep(0.2,4),c(-0.8),c(0.9,-0.4),c(1.2),rep(-0.2,4),rep(0,988)) 
  
  n = 500
  p = 5000
  
  covar <- matrix(rep(0,p*p),nrow=p)
  
  non <-c(1:96,109:4000,4013:5000)
  
  for ( i in non)
    for ( j in non)
      covar[i,j]=0.6^(abs(i-j))
  
  for ( i in 97:100)
    for ( j in 97:100)
      covar[i,j]=0.8
  
  covar[101,101]=0.8
  
  for ( i in 102:104)
    for ( j in 102:104)
      covar[i,j]=0.8
  
  for ( i in 105:108)
    for ( j in 105:108)
      covar[i,j]=0.8
  
  for ( i in 97:100)
    for ( j in 101:101)
      covar[i,j]=0.6
  
  for ( i in 97:100)
    for ( j in 102:104)
      covar[i,j]=0.6
  
  for ( i in 97:100)
    for ( j in 105:108)
      covar[i,j]=0.6
  
  for ( i in 101:101)
    for ( j in 97:100)
      covar[i,j]=0.6
  
  for ( i in 101:101)
    for ( j in 102:104)
      covar[i,j]=0.6
  
  for ( i in 101:101)
    for ( j in 105:108)
      covar[i,j]=0.6
  
  for ( i in 102:104)
    for( j in 97:100)
      covar[i,j]=0.6
  
  for ( i in 102:104)
    for( j in 101:101)
      covar[i,j]=0.6
  
  for ( i in 102:104)
    for( j in 105:108)
      covar[i,j]=0.6
  
  for ( i in 105:108)
    for( j in 97:100)
      covar[i,j]=0.6
  
  for ( i in 105:108)
    for( j in 101:101)
      covar[i,j]=0.6
  
  for ( i in 105:108)
    for( j in 102:104)
      covar[i,j]=0.6
  
  
  for ( i in 4001:4004)
    for ( j in 4001:4004)
      covar[i,j]=0.8
  
  for ( i in 4005:4005)
    for ( j in 4005:4005)
      covar[i,j]=0.8
  
  for ( i in 4006:4007)
    for ( j in 4006:4007)
      covar[i,j]=0.8
  
  for ( i in 4008:4008)
    for ( j in 4008:4008)
      covar[i,j]=0.8
  
  for ( i in 4009:4012)
    for ( j in 4009:4012)
      covar[i,j]=0.8
  
  for ( i in 4001:4004)
    for ( j in 4005:4005)
      covar[i,j]=0.6
  
  for ( i in 4001:4004)
    for ( j in 4006:4007)
      covar[i,j]=0.6
  
  for ( i in 4001:4004)
    for ( j in 4008:4008)
      covar[i,j]=0.6
  
  for ( i in 4001:4004)
    for ( j in 4009:4012)
      covar[i,j]=0.6
  
  for ( i in 4005:4005)
    for ( j in 4001:4004)
      covar[i,j]=0.6
  
  for ( i in 4005:4005)
    for ( j in 4006:4007)
      covar[i,j]=0.6
  
  for ( i in 4005:4005)
    for ( j in 4008:4008)
      covar[i,j]=0.6
  
  for ( i in 4005:4005)
    for ( j in 4009:4012)
      covar[i,j]=0.6
  
  for ( i in 4006:4007)
    for( j in 4001:4004)
      covar[i,j]=0.6
  
  for ( i in 4006:4007)
    for( j in 4005:4005)
      covar[i,j]=0.6
  
  for ( i in 4006:4007)
    for( j in 4008:4008)
      covar[i,j]=0.6
  
  for ( i in 4006:4007)
    for( j in 4009:4012)
      covar[i,j]=0.6
  
  for ( i in 4008:4008)
    for( j in 4001:4004)
      covar[i,j]=0.6
  
  for ( i in 4008:4008)
    for( j in 4005:4005)
      covar[i,j]=0.6
  
  for ( i in 4008:4008)
    for( j in 4006:4007)
      covar[i,j]=0.6
  
  for ( i in 4008:4008)
    for( j in 4009:4012)
      covar[i,j]=0.6
  
  for ( i in 4009:4012)
    for( j in 4001:4004)
      covar[i,j]=0.6
  
  for ( i in 4009:4012)
    for( j in 4005:4005)
      covar[i,j]=0.6
  
  for ( i in 4009:4012)
    for( j in 4006:4007)
      covar[i,j]=0.6
  
  for ( i in 4009:4012)
    for( j in 4008:4008)
      covar[i,j]=0.6
  
  
  for (i in 1:p)
  {
    covar[i,i]=1
  }
  
  system.time(x <- rmultnorm(n,rep(0,p),covar) )
  
  index =seq(1,p)
  fun <- function(j,x)
  {
    t = x[,j]
    n = dim(x)[1]
    d = numeric(n)
    d[t <= -0.44] =0
    d[t>-0.44 & t<0.44]=1
    d[t>=0.44]=2
    list (d=d)
  }
  
  x=matrix(unlist(lapply(index,fun,x)),nrow=n)
  
  res <-rnorm(n = n, 0 , k)
  y <- x%*%beta+res
  snr <- (t(x%*%beta)%*%(x%*%beta)/k^2)^0.5
  list(x=x,y=y,snr=snr)
}

