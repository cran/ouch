weight.matrix <- function (alpha, parsed.tree) {
  N <- parsed.tree$N;
  R <- parsed.tree$R;
  T <- parsed.tree$T;
  ep <- parsed.tree$epochs;
  beta <- parsed.tree$beta;
  W <- matrix(data=0,nrow=N,ncol=R+1);
  W[,1] <- exp(-alpha*T);      
  for (i in 1:N) {
    delta <- diff(exp(alpha*(ep[[i]]-T)));
    for (k in 1:R) {
      W[i,k+1] <- -sum(delta * beta[[i+N*(k-1)]]);
    }
  }
  return(W);
}

