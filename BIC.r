source("sp_scaled.r")

BIC.lambda.w2 <- function (x, y, group, g, n.step,alpha,epsilon) #alpha=1 when lasso
{          
   p <- length(group)
   n <- nrow(x)
   y.c = y - mean(y)
   pos_e = cumsum(group)
   pos_s = pos_e -group + 1

   fun <- function(i,x)
   {
     sigma = t(x[,pos_s[i]:pos_e[i]])%*%x[,pos_s[i]:pos_e[i]]/n
     R = chol(sigma)  #cholesky decomp always exists for X'X/n where dj < n
     x.t = x[,pos_s[i]:pos_e[i]]%*%solve(R)
     return(list(R=R,x.t=x.t))
   }
   result <- lapply(seq(1,p),fun,x) 

   fun <- function(j,result)
   {
     result[[j]]$x.t
   }
   x.t <- matrix(unlist(lapply(seq(1,p),fun,result)),nrow=n)
 
   fun <- function(j,result,y.c)
   {
     t=result[[j]]$x.t
     inner=sum(colSums(t*as.vector(y.c)/n)^2)^0.5/group[j]^0.5
     list(inner=inner)
   }
   lambda.max=max(unlist(lapply(seq(1,p),fun,result,y.c)))
   
   lambda.min <- epsilon*lambda.max
   ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
   lambda1.seq <- numeric(n.step)
   lambda2.seq <- numeric(n.step)
   lh <-numeric(n.step)
   penalty <- numeric(n.step)
   mm <- numeric(n.step)
   bic <- numeric(n.step)
   beta.ini <- rep(0,p)   

   index = seq(1,p-1)
   fun <- function(i,x,pos_s,pos_e)
   {
     a1 <- x[,(pos_s[i]:pos_e[i])]
     a2 <- x[,(pos_s[i+1]:pos_e[i+1])]
     cc <- cancor(a1,a2)
     max(cc$cor)
   }
   w2 <- unlist(lapply(index,fun,x,pos_s,pos_e))   
  
   for ( i in 1: n.step )   
   {
     lambda1.seq[i] <- exp(log(lambda.max)-ss*(i-1))
     lambda2.seq[i] <- lambda1.seq[i]/alpha*(1-alpha)
   }
   
   sp = gslasso_sp(x.t,y.c,w2,group,lambda1.seq,lambda2.seq,T)
   beta = sp$beta

   for ( i in 1: n.step )   
   {
     beta.t <- beta[i,]
     y.h <- x.t%*%as.vector(beta.t)
     lh[i]<- n*log(sum((y.c-y.h)^2)/n)
     mm[i] <- length(beta.t[abs(beta.t)!=0])
     penalty[i] <- mm[i]*(log(n)+2*g*log(p))     #+2*log(d) sum(abs(betam))#
     bic[i] <- lh[i]+penalty[i]
   }
   beta.opt = beta[bic==min(bic),]
   list(bic=bic,lambda.seq=lambda1.seq,lh=lh,penalty=penalty,mm=mm,beta.opt=beta.opt)
}


