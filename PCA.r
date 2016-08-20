PCA <- function(x)
{
  p <- ncol(x)
  R <- cor(x)
  myEig <- eigen(R)
  sum.eigen <- sum(myEig$values)
  csum <- cumsum(myEig$values)/sum.eigen
  index <- seq(1,p)
  d <- min(index[csum > 0.95])
  standardize <- function(x) {(x - mean(x))/sd(x)}
  X <- apply(x, MARGIN=2, FUN=standardize)
  X.t <- X %*% myEig$vectors
  x.s <- X.t[,1:d]
  list(x.s=x.s,d=d)
}

pca.group <- function(i,x,group)
{
  p <- length(group)
  n <- nrow(x)
  pos_e = cumsum(group)
  pos_s = pos_e -group + 1
  pe <- pos_e[i]
  ps <- pos_s[i]

  if ( group[i] != 1 )
  {

    xs = x[,ps:pe]
    pca = PCA(xs)
    x.s = pca$x.s
    d = pca$d
  }
  else if (group[i] == 1)
  {
    x.s = x[,ps:pe]
    d = 1
  }
  list(x.s=x.s,d=d)
}
