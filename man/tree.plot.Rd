\name{tree.plot}
\alias{tree.plot}
\title{Simple plot of phylogenetic tree.}
\description{
  Plot the tree.
}
\usage{
tree.plot(topology, times, names = NULL, regimes = NULL)
}
\arguments{
  \item{topology}{Specification of the topology of the phylogenetic
  tree.  This is in the form of a vector of nonnegative integers, one
  for each node in the tree.  The i-th number is that of the ancestor
  of the i-th node.  The root node is distinguished by the number zero
  (0).}
  \item{times}{A vector of nonnegative numbers, one per node in the
  tree, specifying the time at which each node is located.  The root node
  should be assigned time 0.}
  \item{names}{Optional vector of species names.}
  \item{regimes}{A vector of codes, one for each node in the
  tree, specifying the selective regimes hypothesized to have been
  operative.  Corresponding to each node, enter the code of the regime
  hypothesized for the branch segment terminating in that node.
  For the root node, because it has no branch segment terminating on
  it, the regime specification is irrelevant.}
}
\examples{
data(bimac)
attach(bimac)
tree.plot(ancestor,time,species,OU.LP)
}
\author{Aaron A. King <king at tiem dot utk dot edu>}
\keyword{hplot}
