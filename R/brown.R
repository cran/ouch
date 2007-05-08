brown.fit <- function (data, node, ancestor, times) {
  pt <- parse.tree(node,ancestor,times)
  n <- pt$N
  v <- pt$branch.times
  w <- matrix(data=1,nrow=pt$N,ncol=1)
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  g <- glssoln(w,dat,v)
  theta <- g$coeff
  e <- g$residuals
  sigma <- sqrt((e %*% solve(v,e))/n)
  dim(sigma) <- 1
  u = n * (1 + log(2*pi*sigma*sigma)) + log(det(v))
  dim(u) <- 1
  df <- 2
  list(
       sigma=sigma,
       theta=theta,
       loglik=-u/2,
       deviance=u,
       aic=u+2*df,
       sic=u+log(n)*df,
       df=df
       )
}

brown.dev <- function(n = 1, node, ancestor, times, sigma, theta) {
  pt <- parse.tree(node,ancestor,times)
  v <- pt$branch.times
  x <- rmvnorm(n, rep(theta,dim(v)[1]), as.numeric(sigma^2)*v)
  do.call(
          c,
          apply(
                x,
                1,
                function (z) {
                  y <- rep(NA,length(node))
                  y[pt$term] <- z
                  list(y)
                }
                )
          )
}


