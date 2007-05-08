parse.tree <- function (nodenames, ancestors, times, regime.specs=NULL) {
  nodenames <- as.character(nodenames)
  ancestors <- as.character(ancestors)
  if (!is.valid.ouch.tree(nodenames,ancestors,times,regime.specs))
    stop('the specified tree is not in valid ouch format')
  term <- terminal.twigs(nodenames,ancestors) # get rownumbers of terminal nodes
  N <- length(term)                     # number of terminal nodes
  anc <- ancestor.numbers(nodenames,ancestors)
  bt <- branch.times(anc,times,term)   # absolute times of branch points
  e <- epochs(anc,times,term)
  if (is.null(regime.specs)) {          # useful for BM models
    pt <- list(
               N=N,
               R=0,
               tree.depth=max(times),
               term=term,
               branch.times=bt,
               epochs=e
               )
  } else {                              # useful for Hansen models
    reg <- set.of.regimes(anc,as.factor(regime.specs))
    pt <- list(
               N=N,
               R=length(reg),
               tree.depth = max(times),
               term=term,
               branch.times=bt,
               epochs=e,
               regime.set=reg,
               beta=regimes(anc,times,as.factor(regime.specs),term)
               )
  }
  return(pt)
}

ancestor.numbers <- function (nodenames, ancestors) { # map ancestor names to row numbers
  sapply(ancestors,function(x)charmatch(x,nodenames),USE.NAMES=F)
}

terminal.twigs <- function (nodenames, ancestors) { # numbers of terminal nodes
  which(nodenames %in% setdiff(nodenames,unique(ancestors)))
}

branch.times <- function (ancestors, times, term) {

  N <- length(term)
  tree.depth <- max(times)              # it is assumed that the root node is at time=0

  bt <- matrix(data=0,nrow=N,ncol=N)

  bt[1,1] <- tree.depth
  for (i in 2:N) {
    pedi <- pedigree(ancestors,term[i])
    for (j in 1:(i-1)) {
      pedj <- pedigree(ancestors,term[j])
      for (k in 1:length(pedi)) {
        if (any(pedj == pedi[k])) break
      }
      bt[j,i] <- bt[i,j] <- times[pedi[k]]
    }
    bt[i,i] <- tree.depth
  }
  bt
}

epochs <- function (ancestors, times, term) {
  N <- length(term)
  e <- vector(length=N,mode="list")
  for (k in 1:N) {
    p <- pedigree(ancestors,term[k])
    e[[k]] <- times[p]	
  }
  e
}

set.of.regimes <- function (ancestors, regime.specs) {
  unique(regime.specs[!is.root.node(ancestors)])
}

regimes <- function (ancestors, times, regime.specs, term) {
  N <- length(term)
  reg <- set.of.regimes(ancestors,regime.specs)
  R <- length(reg)
  beta <- vector(R*N, mode="list")
  for (i in 1:N) {
    for (k in 1:R) {
      p <- pedigree(ancestors, term[i])
      n <- length(p)
      beta[[i + N*(k-1)]] <- as.integer(regime.specs[p[1:(n-1)]] == reg[k])
    }
  }    
  beta
}

pedigree <- function (anc, k) {
  p <- k
  k <- anc[k]
  while (!is.root.node(k)) {
    if (k %in% p) stop('this is no tree: circularity detected at node ', k)
    p <- c(p,k)
    k <- anc[k]
  }
  p
}

is.root.node <- function (anc) {
  is.na(anc)
}

