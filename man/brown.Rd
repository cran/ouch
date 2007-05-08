\name{brown.fit}
\alias{brown.fit}
\alias{brown.dev}
\title{Brownian-motion model of evolution along a phylogenetic tree}
\description{
  These functions relate to the Brownian motion model for phylogenetic
  evolution.
  \item{brown}{fits the parameters \eqn{\sigma} and \eqn{\theta} of this
    model to given data.}
  \item{brown.dev}{generates simulated data sets.}
}
\usage{
brown.fit(data, node, ancestor, times)
brown.dev(n = 1, node, ancestor, times, sigma, theta)
}
\arguments{
  \item{data}{Phenotypic data for extant species, i.e., at the
    terminal ends of the phylogenetic tree.}
  \item{node}{Specification of the names of the nodes.}
  \item{ancestor}{Specification of the topology of the phylogenetic
    tree.  This is in the form of a character vector of node names, one
    for each node in the tree.  The i-th name is that of the ancestor
    of the i-th node.  The root node is distinguished by having no
    ancestor (i.e., NA).}
  \item{times}{A vector of nonnegative numbers, one per node in the
    tree, specifying the time at which each node is located.  The root node
    should be assigned time 0.}
  \item{n}{the number of simulated data sets to generate.}
  \item{sigma}{the value of \eqn{\sigma} to be used in the simulations.}
  \item{theta}{the value of \eqn{\theta} to be used in the simulations.}
}
\value{
  \code{brown} returns a list of the following elements:
  \item{sigma}{Maximum likelihood estimate of \eqn{\sigma}.}
  \item{theta}{Maximum likelihood estimate of \eqn{\theta}.}
  \item{loglik}{Log likelihood.}
  \item{deviance}{-2 loglik.}
  \item{aic}{Akaike information criterion.}
  \item{sic}{Schwartz information criterion (=BIC)}
  \item{df}{Number of parameters estimated (= 2).}

  \code{brown.dev} returns a list of \code{n} simulated data sets.
  Each data set corresponds exactly to the data used in the call to
  \code{brown.fit}.
}
\references{
  Butler, M.A. and A.A. King (2004) Phylogenetic comparative analysis: a
  modeling approach for adaptive evolution. American Naturalist, in press.
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
