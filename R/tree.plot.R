tree.plot <- function (topology, times, names = NULL, regimes = NULL) {
  rx <- range(times,na.rm=T);
  rxd <- 0.1*diff(rx);

  if (is.null(regimes))
    regimes <- factor(rep(1,length(topology)));

  levs <- levels(as.factor(regimes));
  palette <- rainbow(length(levs));

  for (r in 1:length(levs)) {
    y <- tree.layout(topology);
    x <- times;
    f <- which(topology > 0 & regimes == levs[r]);
    pp <- topology[f];
    X <- array(data=c(x[f], x[pp], rep(NA,length(f))),dim=c(length(f),3));
    Y <- array(data=c(y[f], y[pp], rep(NA,length(f))),dim=c(length(f),3));
    oz <- array(data=1,dim=c(2,1));
    X <- kronecker(t(X),oz);
    Y <- kronecker(t(Y),oz);
    X <- X[2:length(X)];
    Y <- Y[1:(length(Y)-1)];
    C <- rep(palette[r],length(X));
    if (r > 1) par(new=T);
    par(yaxt='n')
    plot(X,Y,type='l',col=C,xlab='time',ylab='',xlim = rx + c(-rxd,rxd),ylim=c(0,1));
    if (!is.null(names))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],names[f],pos=4);
  }
}

tree.layout <- function (topology) {
  root <- which(topology==0);
  return(arrange.tree(root,topology));
}

arrange.tree <- function (root, topology) {
  k <- which(topology==root);
  n <- length(k);
  reltree <- rep(0,length(topology));
  reltree[root] <- 0.5;
  p <- list(NULL);
  if (n > 0) {
    m <- rep(0,n);
    for (j in 1:n) {
      p[[j]] <- arrange.tree(k[j],topology);
      m[j] <- length(which(p[[j]] != 0));
    }
    cm <- c(0,cumsum(m));
    for (j in 1:n) {
      reltree <- reltree + (cm[j]/sum(m))*(p[[j]] != 0) + (m[j]/sum(m))*p[[j]];
    }
  }
  return(reltree);
}
