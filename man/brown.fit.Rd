\name{brown.fit}
\alias{brown.fit}
\alias{brown.dev}
\title{Brownian-motion model of evolution along a phylogenetic tree}
\description{
  These functions relate to the Brownian motion model for phylogenetic
  evolution.
  \item{brown}{fits the parameters sigma and theta of this model to
    given data.}
  \item{brown.dev}{generates simulated data sets.}
}
\usage{
brown.fit(data, topology, times)
brown.dev(n = 1, topology, times, sigma, theta)
}
\arguments{
  \item{data}{Phenotypic data for extant species, i.e., at the
    terminal ends of the phylogenetic tree.}
  \item{topology}{Specification of the topology of the phylogenetic
    tree.  This is in the form of a vector of nonnegative integers, one
    for each node in the tree.  The i-th number is that of the ancestor
    of the i-th node.  The root node is distinguished by the number zero
    (0).}
  \item{times}{A vector of nonnegative numbers, one per node in the
    tree, specifying the time at which each node is located.  The root node
    should be assigned time 0.}
  \item{n}{the number of simulated data sets to generate.}
  \item{sigma}{the value of sigma to be used in the simulations.}
  \item{theta}{the value of theta to be used in the simulations.}
}
\value{
  brown returns a list of the following elements:
  \item{sigma}{Maximum likelihood estimate of sigma.}
  \item{theta}{Maximum likelihood estimate of theta.}
  \item{u}{-2 log likelihood.}
  \item{aic}{Akaike information criterion.}
  \item{sic}{Schwartz information criterion (=BIC)}
  \item{df}{Number of parameters estimated (= 2).}

  brown.dev returns a data frame containing simulated data sets.  Each
  realization is a row.
}
\references{
  Butler, M.A. and A.A. King (2004) Phylogenetic comparative analysis: a
  modeling approach for adaptive evolution. American Naturalist, in press.
}
\author{Aaron A. King <king at tiem dot utk dot edu>}
\keyword{models}
