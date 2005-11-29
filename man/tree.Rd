\name{tree.plot}
\alias{tree.plot}
\alias{is.valid.ouch.tree}
\title{Validating and plotting phylogenetic trees.}
\description{
  Validate or plot the phylogenetic tree in ouch format.
}
\usage{
is.valid.ouch.tree(node, ancestor, times, regimes = NULL)
tree.plot(node, ancestor, times, names = NULL, regimes = NULL)
}
\arguments{
  \item{node}{A character vector giving the name of each node.}
  \item{ancestor}{Specification of the topology of the phylogenetic
    tree.  This is in the form of a character vector naming the immediate
    ancestor of each node.  In particular, the i-th name is that of the
    ancestor of the i-th node.  The root node is distinguished by having
    no ancestor (i.e., NA).}
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
\details{
  \code{tree.plot} makes a simple plot of the phylogenetic tree.  Labels
  and/or a selective-regime-based coloring scheme are optional.

  \code{is.valid.ouch.tree} performs several checks to make sure that
  the tree is valid.  These include checks to ensure that all the vector
  arguments are of the same length, that the nodes have unique names,
  that there is a unique root, that there is at least one terminal node,
  that every node's ancestor is in fact part of the tree, that there are
  no cycles in the 'tree', and that the tree is connected.  It returns
  TRUE if the tree is valid, FALSE otherwise and gives diagnostic
  warnings.
}
\examples{
data(bimac)
attach(bimac)
is.valid.ouch.tree(node,ancestor,time,OU.LP)
tree.plot(node,ancestor,time,species,OU.LP)
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{hplot}
