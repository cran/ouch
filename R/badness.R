badness <- function (alpha, data, parsed.tree) {
  a <- exp(alpha);
  n <- length(data);
  w <- weight.matrix(a, parsed.tree);
  v <- scaled.covariance.matrix(a, parsed.tree);
  g <- glssoln(w,data,v);
  e <- g$residuals;
  sigmasq <- (e %*% solve(v,e)) / n;
  dim(sigmasq) <- 1;
  u <- n * (1 + log(2*pi*sigmasq)) + log(det(v));
  dim(u) <- 1;
  return(u);
}

