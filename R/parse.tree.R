parse.tree <- function (topology, times, regime.specs=NULL) {
  if ((length(topology) != length(times)) 
       || ((!is.null(regime.specs))
          && (length(topology) != length(regime.specs)))
      ) {
    warning('invalid tree')
    return(NULL);
  }
  term <- terminal.twigs(topology,times);
  N <- length(term);
  T <- max(times);
  bt <- branch.times(topology,times);
  e <- epochs(topology,times,term);
  if (is.null(regime.specs)) {
    pt <- list(N=N,T=T,term=term,branch.times=bt,epochs=e);
  } else {
    reg <- set.of.regimes(topology,regime.specs);    
    R <- length(reg);
    b <- regimes(topology, times, regime.specs, term);
    pt <- list(N=N,R=R,T=T,term=term,branch.times=bt,epochs=e,regime.set=reg,beta=b);
  }
  return(pt);
}

