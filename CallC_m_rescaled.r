dyn.load("gslasso.m.so")

standard <- function (x) {
  p <- ncol(x)
  n <- nrow(x)
  x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
  x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}


gslasso <- function(x,y,w2,G,lambda1,lambda2,orthonormalize=F)
{
  n <- length(y)
  totalp <- dim(x)[2]
  p <- length(G)
  iter <- 200
  
  pos_e <- cumsum(G)
  pos_s <- pos_e - G + 1
  
  if ( orthonormalize == F)
  {
    fun <- function(i,x)
    {
      sigma = t(x[,pos_s[i]:pos_e[i]])%*%x[,pos_s[i]:pos_e[i]]/n
      R = chol(sigma) 
      x.t = x[,pos_s[i]:pos_e[i]]%*%solve(R)
      return(list(R=R,x.t=x.t))
    }
    result <- lapply(seq(1,p),fun,x) 
    
    fun <- function(j,result)
    {
      result[[j]]$x.t
    }
    x.t <- matrix(unlist(lapply(seq(1,p),fun,result)),nrow=n)
  }
  else if ( orthonormalize == T)
    x.t = x 
  
  beta <- numeric(totalp)
  param <- c(n, p, totalp,iter,max(G)) 
  epsilon <- 1E-10
  fit <- .C("Gslasso", y=as.double(y),x=as.double(t(x.t)),G=as.integer(G),param=as.integer(param),lambda=as.double(c(lambda1,lambda2)),w2=as.double(w2),epsilon=as.double(epsilon),beta=as.double(beta))
  
  beta <- numeric()
  b = fit$beta
  
  index <- seq(1,p)
  fun <- function(i,w2,G)
  {
    if ( i == 1 )
      a = w2[i]*max(G)/G[i]
    else if ( i == p )
      a = w2[i-1]*max(G)/G[i]
    else 
      a = (w2[i-1]+w2[i])*max(G)/G[i]
  }
  c <- unlist(lapply(index,fun,w2,G))
  c <- rep(c,G)
  
  for ( j in 1:p )
  {
    if ( orthonormalize == F)
    {
      beta <- c(beta,solve(result[[j]]$R,b[pos_s[j]:pos_e[j]]))
    }
    else if (orthonormalize == T)
    {
      beta <- b
    }
  }
  
  beta.resc <- beta*(1+lambda2*c)
  
  list(beta=beta,beta.resc=beta.resc,b=b)
}
