tree.plot <- function (node, ancestor, times, names = NULL, regimes = NULL) {

  node <- as.character(node)
  ancestor <- as.character(ancestor)
  if (!is.valid.ouch.tree(node,ancestor,times,regimes))
    stop("the tree is not in valid ouch format");

  rx <- range(times,na.rm=T)
  rxd <- 0.1*diff(rx)

  anc <- ancestor.numbers(node,ancestor)

  if (is.null(regimes))
    regimes <- factor(rep(1,length(anc)))

  levs <- levels(as.factor(regimes))
  palette <- rainbow(length(levs))

  for (r in 1:length(levs)) {
    y <- tree.layout(anc)
    x <- times
    f <- which(!is.root.node(anc) & regimes == levs[r])
    pp <- anc[f]
    X <- array(data=c(x[f], x[pp], rep(NA,length(f))),dim=c(length(f),3))
    Y <- array(data=c(y[f], y[pp], rep(NA,length(f))),dim=c(length(f),3))
    oz <- array(data=1,dim=c(2,1))
    X <- kronecker(t(X),oz)
    Y <- kronecker(t(Y),oz)
    X <- X[2:length(X)]
    Y <- Y[1:(length(Y)-1)]
    C <- rep(palette[r],length(X))
    if (r > 1) par(new=T)
    par(yaxt='n')
    plot(X,Y,type='l',col=C,xlab='time',ylab='',xlim = rx + c(-rxd,rxd),ylim=c(0,1))
    if (!is.null(names))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],names[f],pos=4)
  }
}

tree.layout <- function (anc) {
  root <- which(is.root.node(anc))
  arrange.tree(root,anc)
}

arrange.tree <- function (root, anc) {
  k <- which(anc==root)
  n <- length(k)
  reltree <- rep(0,length(anc))
  reltree[root] <- 0.5
  p <- list()
  if (n > 0) {
    m <- rep(0,n)
    for (j in 1:n) {
      p[[j]] <- arrange.tree(k[j],anc)
      m[j] <- length(which(p[[j]] != 0))
    }
    cm <- c(0,cumsum(m))
    for (j in 1:n) {
      reltree <- reltree + (cm[j]/sum(m))*(p[[j]] != 0) + (m[j]/sum(m))*p[[j]]
    }
  }
  reltree
}


is.valid.ouch.tree <- function (node, ancestor, times, regimes = NULL) {
  valid <- TRUE
  node <- as.character(node)
  ancestor <- as.character(ancestor)
  n <- length(node)
  if (length(unique(node)) != n) {
    warning('node names must be unique')
    valid <- FALSE
  }
  if (
      (length(ancestor) != n) ||
      (length(times) != n)
      ) {
    warning('invalid tree: all columns must be of the same length')
    valid <- FALSE
  }
  if (!is.null(regimes) && (length(regimes) != n)) {
    warning('regimes must be the same length as the other columns')
    valid <- FALSE
  }
  root <- which(is.root.node(ancestor))
  if (length(root) != 1) {
    warning('the tree must have a unique root node, designated by its having ancestor = NA')
    return(FALSE)
  }
  term <- as.list(terminal.twigs(node,ancestor))
  if (length(term) <= 0) {
    warning("there ought to be at least one terminal node, don't you think?")
    valid <- FALSE
  }
  outs <- which((!is.root.node(ancestor) & !(ancestor %in% node)))
  if (length(outs) > 0) {
    str <- sprintf("the ancestor of node %s is not in the tree\n", node[outs])
    warning(str,call.=F)
    valid <- FALSE
  }
  anc <- ancestor.numbers(node,ancestor)
  ck <- all(
            sapply(
                   1:n,
                   function(x) {
                     good <- root %in% pedigree(anc,x)
                     if (!good) {
                       str <- sprintf("node %s is disconnected", node[x])
                       warning(str,call.=F)
                     }
                     good
                   }
                   )
            )
  valid && ck
}

