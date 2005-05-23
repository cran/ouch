\name{hansen.fit}
\alias{hansen.fit}
\alias{hansen.dev}
\title{Hansen model of evolution along a phylogenetic tree}
\description{
  These functions relate to the Hansen model for phylogenetic evolution.
  \item{hansen.fit}{fits the parameters alpha, sigma, and theta of this
    model.}
  \item{hansen.dev}{generates random deviates from the Hansen model.}
}
\usage{
hansen.fit(data, node, ancestor, times, regimes,
           guess = 0, interval = c(0.001, 20), tol = 1e-12)
hansen.dev(n = 1, node, ancestor, times, regimes, alpha, sigma, theta)
}
\arguments{
  \item{data}{Phenotypic data for extant species, i.e., at the
    terminal ends of the phylogenetic tree.}
  \item{node}{Character vector of names of the nodes.}
  \item{ancestor}{Specification of the topology of the phylogenetic
    tree.  This is in the form of a character vector of node names, one
    for each node in the tree.  The i-th name is that of the ancestor of
    the i-th node.  The root node is distinguished by having no ancestor
    (i.e., NA).}
  \item{times}{A vector of nonnegative numbers, one per node in the
    tree, specifying the time at which each node is located.  The root node
    should be assigned time 0.}
  \item{regimes}{A vector of codes, one for each node in the
    tree, specifying the selective regimes hypothesized to have been
    operative.  Corresponding to each node, enter the code of the regime
    hypothesized for the branch segment terminating in that node.
    For the root node, because it has no branch segment terminating on
    it, the regime specification is irrelevant.}
  \item{guess}{(Optional).  An initial guess for the selection
    strength parameter alpha.  This is used to initialize the
    optimization algorithm.}
  \item{interval}{The interval which will be searched for the optimal
    alpha.  By default, 0.001 < alpha < 20.}
  \item{tol}{Convergence tolerance.}
  \item{n}{Number of pseudorandom data sets to generate.}
  \item{alpha}{Value of alpha to use.}
  \item{sigma}{Value of sigma to use.}
  \item{theta}{Value of theta to use.}
}
\value{
  hansen.fit returns a list containing the following elements:
  \item{alpha}{Maximum likelihood estimate of alpha.  Note that if
    alpha lies against one of the constraints (see \code{interval}
    above), then this is not a maximum-likelihood estimate.} 
  \item{sigma}{Maximum likelihood estimate of sigma.}
  \item{theta}{Maximum likelihood estimate of theta.}
  \item{u}{-2 log likelihood.}
  \item{aic}{Akaike information criterion.}
  \item{sic}{Schwartz information criterion (=BIC)}
  \item{df}{Number of parameters estimated (= 3 + number of regimes).}

  hansen.dev returns a data frame containing simulated data sets.  Each
  realization is a row.
}
\references{
  Butler, M.A. and A.A. King (2004) Phylogenetic comparative analysis: a
  modeling approach for adaptive evolution. American Naturalist, in press.
}
\author{Aaron A. King <king at tiem dot utk dot edu>}
\keyword{models}
