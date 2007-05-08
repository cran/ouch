rmvnorm <- function (n = 1, mu, sigma) {
  p <- length(mu)
  if (!all(dim(sigma) == c(p,p)))
    stop("incompatible arguments")
  cf <- chol(sigma,pivot=F)
  X <- matrix(mu,n,p,byrow=T)+matrix(rnorm(p*n),n)%*%cf
  if (n == 1) {
    drop(X)
  }  else {
    X
  }
}
