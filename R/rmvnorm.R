rmvnorm <- function (n = 1, mu, sigma, tol = 1e-06) {
  p <- length(mu);
  if (!all(dim(sigma) == c(p,p)))
    stop("incompatible arguments");
  cf <- chol(sigma,pivot=F);
  X <- matrix(mu,n,p,byrow=T) + matrix(rnorm(p*n),n) %*% cf;
  if (n == 1) {
    return(drop(X));
  }  else {
    return(X);
  }
}

