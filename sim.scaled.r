source("BIC.r")
source("data.r")
source("PCA.r")

match.fun <- function(betaf)
{  
  a=c(97:105,4001:4009)  
  b=seq(1:5000)
  b=b[betaf!=0]
  matched=match(a,b,nomatch=-1)
  tp <- length(matched[matched!=-1])
  fn <- 18 -tp
  fp <- length(b) - tp
  tn <- 5000 - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

match.fun.group <- function(betaf,group)
{  
  pos_e = cumsum(group)
  pos_s = pos_e -group + 1
  p = length(group)  
  tp = 0
  for ( i in c(25:28,1002:1006) )
  {
    if ( betaf[pos_s[i]]!=0 )
      tp = tp +1
  }

  fp = 0
  for ( i in c(1:24,29:1001,1007:1253) )
  {
    if ( betaf[pos_s[i]]!=0 )
      fp = fp +1
  }

  fn <- 9 -tp
  tn <- 1253 - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

sim <- function(b,alpha,n.step,epsilon,g)
{
  tp <- numeric()
  tp1 <- numeric()
  tp2 <- numeric()
  tn <- numeric() 
  fp <- numeric()
  fn <- numeric()
  lambda.opt.seq <- numeric()
  lambda.min.seq <- numeric() 
  lambda.max.seq <- numeric()
  group <- c(rep(4,24),4,1,3,4,rep(4,973),4,1,2,1,4,rep(4,247))

  for ( i in 1:b )
  {
    data=generate(1.5) 
    x = data$x
    y = data$y   

    n = dim(x)[1]    

    index=seq(1,length(group))
    tmp=lapply(index,pca.group,x,group)

    fun <- function(i,tmp)
    {
      x = tmp[[i]]$x.s
    }
    xnew=matrix(unlist(lapply(index,fun,tmp)),nrow=n)

    fun <- function(i,tmp)
    {
     x = tmp[[i]]$d
    }
    group.new=unlist(lapply(index,fun,tmp))

    opt=BIC.lambda.w2(xnew,y,group.new,1,g,n.step,alpha,epsilon)
    lambda.seq = opt$lambda.seq
    lambda.opt = lambda.seq[opt$bic==min(opt$bic)]
    lambda.opt.seq <- c(lambda.opt.seq,lambda.opt)
    lambda.min.seq <- c(lambda.min.seq,min(lambda.seq))
    lambda.max.seq <- c(lambda.max.seq,max(lambda.seq))
    beta.opt=opt$beta.opt
    a=match.fun.group(beta.opt,group.new)
    tp = c(tp,a$tp)
    tn = c(tn,a$tn)
    fp = c(fp,a$fp)
    fn = c(fn,a$fn)
  }
  list(tp=tp,tn=tn,fp=fp,fn=fn,lambda.opt.seq=lambda.opt.seq,lambda.min.seq=lambda.min.seq,lambda.max.seq=lambda.max.seq)
}