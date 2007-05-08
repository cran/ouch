\name{hansen.fit}
\alias{hansen.fit}
\alias{hansen.dev}
\alias{hansen.prof}
\title{Hansen model of evolution along a phylogenetic tree}
\description{
  These functions relate to the Hansen model for phylogenetic evolution.
  \item{hansen.fit}{fits the parameters \eqn{\alpha}, \eqn{\sigma}, and
    \eqn{\theta}  of this model.}
  \item{hansen.prof}{computes the profile likelihood of the \eqn{\alpha}
    parameter.}
  \item{hansen.dev}{generates random deviates from the Hansen model.}
}
\usage{
hansen.fit(data, node, ancestor, times, regimes = NULL,
           interval = c(0, 100), tol = 1e-12)
hansen.prof(alpha, data, node, ancestor, times, regimes = NULL)
hansen.dev(n = 1, node, ancestor, times, regimes = NULL,
           alpha, sigma, theta)
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
  \item{interval}{The interval which will be searched for the optimal
    \eqn{\alpha}.  By default, \eqn{0.001 < \alpha < 20}.}
  \item{tol}{Convergence tolerance.}
  \item{n}{Number of pseudorandom data sets to generate.}
  \item{alpha}{Value of \eqn{\alpha} to use.}
  \item{sigma}{Value of \eqn{\sigma} to use.}
  \item{theta}{Value of \eqn{\theta} to use.}
}
\value{
  \code{hansen.fit} returns a list containing the following elements:
  \item{alpha}{Maximum likelihood estimate of \eqn{\alpha}.  Note that
    if \eqn{\alpha} lies against one of the constraints (see
    \code{interval} above), then this is not a maximum-likelihood
    estimate.}
  \item{sigma}{Maximum likelihood estimate of \eqn{\sigma}.}
  \item{theta}{Maximum likelihood estimate of \eqn{\theta}.}
  \item{loglik}{Log likelihood.}
  \item{deviance}{-2 loglik.}
  \item{aic}{Akaike information criterion.}
  \item{sic}{Schwartz information criterion (=BIC)}
  \item{df}{Number of parameters estimated (= 3 + number of regimes).}

  \code{hansen.prof} allows one to compute the profile log likelihood
  across a range of \eqn{\alpha}.  Specifically, at the values of
  \eqn{\alpha} specified in \code{alpha}, the log likelihood of the given
  data under the given regimes is computed.  \code{hansen.prof} returns
  a named list containing the following elements:
  \item{alpha}{the specified \code{alpha}}
  \item{loglik}{the log likelihood}
  \item{deviance}{the deviance (\eqn{-2 \log L})}
  \item{aic}{the Akaike information criterion value}
  \item{sic}{the Schwartz information criterion value}
  Note that when \eqn{\alpha=0} exactly, the computed log likelihood
  does not agree with the Brownian motion model.

  \code{hansen.dev} returns a list of \code{n} simulated data sets.
  Each data set corresponds exactly to the data used in the call to
  \code{hansen.fit}.
}
\details{
  The vector \code{regimes} should be of class factor.  If
  \code{regimes} is unspecified or NULL, all lineages are assumed to be
  evolving under a single, global OU process with a global optimum.
  In this case, rather than estimate the character state at the root
  node, the algorithm assumes that the character state at the root value
  follows the stationary distribution for the OU process.  In general,
  it is impossible to identify both the root character state and the
  global optimum using contemporaneous data.
}
\references{
  Butler, M.A. and A.A. King (2004) Phylogenetic comparative analysis: a
  modeling approach for adaptive evolution. American Naturalist
  164:683-695, 2004.
}
\author{Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
