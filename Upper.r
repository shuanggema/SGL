source("CallC_m_rescaled.r")

group.id <- function(beta,group)
{
  pos_e = cumsum(group)
  pos_s = pos_e -group + 1
  group.s <- numeric()
  for ( i in 1:length(group) )
  {
    if ( beta[pos_s[i]] != 0 )
      group.s <- c(group.s, i)
  }
  list(group.s=group.s)
}


upper <- function (x, y, group, pre.s, alpha,epsilon,max) #alpha=1 when lasso
{          
   p <- length(group)
   n <- nrow(x)
   y.c = y - mean(y)
   pos_e = cumsum(group)
   pos_s = pos_e -group + 1

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
 
   fun <- function(j,result,y.c)
   {
     t=result[[j]]$x.t
     inner=sum(colSums(t*as.vector(y.c)/n)^2)^0.5/group[j]^0.5
     list(inner=inner)
   }
   lambda.max=max(unlist(lapply(seq(1,p),fun,result,y.c)))

   index = seq(1,p-1)
   fun <- function(i,x,pos_s,pos_e)
   {
     a1 <- x[,(pos_s[i]:pos_e[i])]
     a2 <- x[,(pos_s[i+1]:pos_e[i+1])]
     cc <- cancor(a1,a2)
     max(cc$cor)
   }
   w2 <- unlist(lapply(index,fun,x,pos_s,pos_e))     


   upper.l <- log(lambda.max)
   lower.l <- log(lambda.max*epsilon)
   lambda1 = exp((lower.l + upper.l)/2)
   lambda2 = lambda1/alpha*(1-alpha)
   sp = gslasso(x.t,y.c,w2,group,lambda1,lambda2,T)
   beta = sp$beta   
   l = 0 

   repeat{ 
      group.s = group.id(beta,group)$group.s
      if ( length(group.s) == pre.s | l == max) 
      {
         break
      }  
      else if ( length(group.s) > pre.s )
      {
         lower.l <- log(lambda1)
         lambda1 <- exp((lower.l + upper.l)/2)
         lambda2 <- lambda1/alpha*(1-alpha)
         sp = gslasso(x.t,y.c,w2,group,lambda1,lambda2,T)
         beta = sp$beta
      }
      else if ( length(group.s) < pre.s )
      {
         upper.l <- log(lambda1)
         lambda1 <- exp((lower.l + upper.l)/2)
         lambda2 <- lambda1/alpha*(1-alpha)
         sp = gslasso(x.t,y.c,w2,group,lambda1,lambda2,T)
         beta = sp$beta
      }
      l = l + 1
   }

   beta.b <- numeric()
   for ( j in 1:p )
   {
      beta.b <- c(beta.b,solve(result[[j]]$R,beta[pos_s[j]:pos_e[j]]))
   }
   
   list(lambda1=lambda1,lambda2=lambda2,group.s=group.s,beta=beta.b,beta.o=beta,w2=w2,lambda.max=lambda.max,x.t=x.t)
}

