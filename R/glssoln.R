glssoln <- function(a, x, v, tol = 1e-12) {
  n <- length(x);
  vh <- t(chol(v));
  s <- svd(forwardsolve(vh,a));
                                        #   Can we be certain that the singular values are sorted in
                                        #   decreasing order?  (Probably not)
                                        #   k <- order(s$d,decreasing=T);
  svals <- s$d[s$d > tol * max(s$d)];
  r <- length(svals);
  svals <-  diag(1/svals,nrow=r,ncol=r);
  y <- (s$v[,1:r] %*% (svals %*% t(s$u[,1:r]))) %*% forwardsolve(vh,x);
  e <- a %*% y - x;
  dim(y) <- dim(y)[1];
  dim(e) <- n;
  return(list(coeff=y,residuals=e));
}

