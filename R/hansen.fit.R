hansen.fit <- function (data, topology, times, regimes, guess=0, interval=c(0.001,20), tol=1e-12) {
  pt <- parse.tree(topology,times,regimes);
  n <- pt$N;
  dat <- data[!is.na(data)];
  r <- optimize(badness,interval=log(interval),
                lower=log(interval[1]),upper=log(interval[2]),
                tol=tol,maximum=F,dat,pt);
  alpha = exp(r$minimum);

  w <- weight.matrix(alpha, pt);
  v <- scaled.covariance.matrix(alpha, pt);
  g <- glssoln(w,dat,v);
  theta <- g$coeff;
  e <- g$residuals;
  sigma <- sqrt((e %*% solve(v,e))/n);
  dim(sigma) <- 1;
  u = r$objective;
  dim(u) <- 1;
  df <- pt$R+3;
  return(list(alpha=alpha,sigma=sigma,theta=theta,u=u,aic=u+2*df,sic=u+log(n)*df,df=df));
}

