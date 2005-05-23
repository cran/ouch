hansen.fit <- function (data, node, ancestor, times, regimes,
                        guess=0, interval=c(0.001,20), tol=1e-12) {
  pt <- parse.tree(node,ancestor,times,regimes);
  n <- pt$N;
  dat <- data[pt$term]
  no.dats <- which(is.na(dat))
  if (length(no.dats) > 0)
    stop("Missing data on terminal nodes: ",node[pt$term[no.dats]])
  r <- optimize(badness,interval=log(interval),
                lower=log(interval[1]),upper=log(interval[2]),
                tol=tol,maximum=F,dat,pt);
  alpha = exp(r$minimum);

  w <- weight.matrix(alpha, pt);
  v <- scaled.covariance.matrix(alpha, pt);
  g <- glssoln(w,dat,v);
  theta <- g$coeff;
  names(theta) <- paste('theta',c('0',as.character(pt$regime.set)),sep='.')
  e <- g$residuals;
  sigma <- sqrt((e %*% solve(v,e))/n);
  dim(sigma) <- 1;
  u = r$objective;
  dim(u) <- 1;
  df <- pt$R+3;
  return(list(alpha=alpha,sigma=sigma,theta=theta,u=u,aic=u+2*df,sic=u+log(n)*df,df=df));
}

hansen.dev <- function(n = 1, node, ancestor, times, regimes, alpha, sigma, theta) {
  pt <- parse.tree(node,ancestor,times,regimes);
  w <- weight.matrix(alpha, pt);
  v <- scaled.covariance.matrix(alpha, pt);
  x <- rmvnorm(n, as.vector(w %*% theta), as.numeric(sigma^2)*v);
  return(data.frame(x))
}

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

weight.matrix <- function (alpha, parsed.tree) {
  N <- parsed.tree$N;
  R <- parsed.tree$R;
  tree.depth <- parsed.tree$tree.depth;
  ep <- parsed.tree$epochs;
  beta <- parsed.tree$beta;
  W <- matrix(data=0,nrow=N,ncol=R+1);
  W[,1] <- exp(-alpha*tree.depth);      
  for (i in 1:N) {
    delta <- diff(exp(alpha*(ep[[i]]-tree.depth)));
    for (k in 1:R) {
      W[i,k+1] <- -sum(delta * beta[[i+N*(k-1)]]);
    }
  }
  return(W);
}

scaled.covariance.matrix <- function (alpha, parsed.tree) {
  tree.depth <- parsed.tree$tree.depth;
  bt <- parsed.tree$branch.times;			 
  if (alpha == 0) {
    V <- bt;
  } else {
    a <- 2*alpha;
    V <- exp(-a*tree.depth) * expm1(a*bt) / a;
  }
}

