brown.fit <- function (data, topology, times) {
  pt <- parse.tree(topology,times);
  n <- pt$N;
  v <- pt$branch.times;
  w <- matrix(data=1,nrow=pt$N,ncol=1);
  dat <- data[!is.na(data)];
  g <- glssoln(w,dat,v);
  theta <- g$coeff;
  e <- g$residuals;
  sigma <- sqrt((e %*% solve(v,e))/n);
  dim(sigma) <- 1;
  u = n * (1 + log(2*pi*sigma*sigma)) + log(det(v));
  dim(u) <- 1;
  df <- 2;
  return(list(sigma=sigma,theta=theta,u=u,aic=u+2*df,sic=u+log(n)*df,df=df));
}

