hansen.dev <- function(n = 1, topology, times, regimes, alpha, sigma, theta) {
  pt <- parse.tree(topology,times,regimes);
  w <- weight.matrix(alpha, pt);
  v <- scaled.covariance.matrix(alpha, pt);
  x <- rmvnorm(n, as.vector(w %*% theta), as.numeric(sigma^2)*v);
  return(data.frame(x))
}
