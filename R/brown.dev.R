brown.dev <- function(n = 1, topology, times, sigma, theta) {
  pt <- parse.tree(topology,times);
  v <- pt$branch.times;
  x <- rmvnorm(n, rep(theta,dim(v)[1]), as.numeric(sigma^2)*v);
  return(data.frame(x))
}


