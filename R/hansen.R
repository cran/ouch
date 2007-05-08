hansen.fit <- function (data, node, ancestor, times,
                        regimes = NULL,
                        interval = c(0,100),
                        tol = 1e-12) {
  pt <- parse.tree(node,ancestor,times,regimes)
  n <- pt$N
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  r <- optimize(
                f=badness,
##                interval=log(interval),
                interval=interval,
                tol=tol,
                maximum=F,
                data=dat,
                parsed.tree=pt
                )
##  alpha <- exp(r$minimum)
  alpha <- r$minimum
  w <- weight.matrix(alpha,pt)
  v <- scaled.covariance.matrix(alpha,pt)
  g <- glssoln(w,dat,v)
  theta <- g$coeff
  if (pt$R>0) {
    names(theta) <- c('0',as.character(pt$regime.set))
  }
  e <- g$residuals
  sigma <- sqrt((e %*% solve(v,e))/n)
  dim(sigma) <- 1
  u = r$objective
  dim(u) <- 1
  df <- pt$R+3
  list(
       alpha=alpha,
       sigma=sigma,
       theta=theta,
       loglik=-u/2,
       deviance=u,
       aic=u+2*df,
       sic=u+log(n)*df,
       df=df
       )
}

hansen.dev <- function(n = 1, node, ancestor, times, regimes = NULL, alpha, sigma, theta) {
  pt <- parse.tree(node,ancestor,times,regimes)
  w <- weight.matrix(alpha,pt)
  v <- scaled.covariance.matrix(alpha,pt)
  x <- rmvnorm(n,as.vector(w%*%theta),as.numeric(sigma^2)*v)
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

hansen.prof <- function (alpha,
                         data, node, ancestor, times,
                         regimes = NULL) {
  pt <- parse.tree(node,ancestor,times,regimes)
  n <- pt$N
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  u <- sapply(alpha,badness,data=dat,parsed.tree=pt)
  df <- pt$R+3
  list(
       alpha=alpha,
       loglik=-u/2,
       deviance=u,
       aic=u+2*df,
       sic=u+log(n)*df
       )
}

badness <- function (alpha, data, parsed.tree) {
##  a <- exp(alpha)
  a <- alpha
  n <- length(data)
  w <- weight.matrix(a,parsed.tree)
  v <- scaled.covariance.matrix(a,parsed.tree)
  g <- glssoln(w,data,v)
  e <- g$residuals
  sigmasq <- (e %*% solve(v,e)) / n
  dim(sigmasq) <- 1
  u <- n * (1 + log(2*pi*sigmasq)) + log(det(v))
  dim(u) <- 1
  u
}

weight.matrix <- function (alpha, parsed.tree) {
  N <- parsed.tree$N
  R <- parsed.tree$R
  if (R > 0) {
    tree.depth <- parsed.tree$tree.depth
    ep <- parsed.tree$epochs
    beta <- parsed.tree$beta
    w <- matrix(data=0,nrow=N,ncol=R+1)
    w[,1] <- exp(-alpha*tree.depth)
    for (i in 1:N) {
      delta <- diff(exp(alpha*(ep[[i]]-tree.depth)))
      for (k in 1:R) {
        w[i,k+1] <- -sum(delta * beta[[i+N*(k-1)]])
      }
    }
  } else {
    w <- matrix(data=1,nrow=parsed.tree$N,ncol=1)
  }
  w
}

scaled.covariance.matrix <- function (alpha, parsed.tree) {
  tree.depth <- parsed.tree$tree.depth
  bt <- parsed.tree$branch.times
  if (alpha == 0) {
    v <- bt
  } else {
    a <- 2*alpha
    if (parsed.tree$R > 0) {
      v <- exp(-a*tree.depth)*expm1(a*bt)/a
    } else {
      v <- exp(-a*(tree.depth-bt))/a
    }
  }
  v
}

