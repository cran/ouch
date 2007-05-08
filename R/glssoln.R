glssoln <- function (a, x, v, tol = sqrt(.Machine$double.eps)) {
  n <- length(x)
  vh <- t(chol(v))
  s <- svd(forwardsolve(vh,a))
  inds <- s$d > tol*max(s$d)
  svals <- s$d[inds,drop=FALSE]
  r <- length(svals)
  svals <-  diag(1/svals,nrow=r,ncol=r)
  y <- (s$v[,inds,drop=FALSE] %*% (svals %*% t(s$u[,inds,drop=FALSE]))) %*% forwardsolve(vh,x)
  e <- a %*% y - x
  dim(y) <- dim(y)[1]
  dim(e) <- n
  list(
       coeff=y,
       residuals=e
       )
}
