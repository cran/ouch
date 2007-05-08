\name{bimac}
\alias{bimac}
\docType{data}
\title{Anolis bimaculatus lizard size data.}
\description{
  This is the Anolis bimaculatus dataset used in Butler & King
  2004. It is used to test a hypothesis of character displacement
  using an interspecific dataset of body sizes and current data on
  sympatry/allopatry.  The data frame consists of the following
  columns: \code{species} which are species names, \code{size} which
  is the phenotypic data, and the variables \code{ancestor} and
  \code{time} which specify the topology of the phylogeny and the
  location of the nodes in time, respectively. The columns
  \code{OU.1}, \code{OU.3}, \code{OU.4}, and \code{OU.LP} specify the
  hypothesized selective regimes (there are four alternatives in this
  example).  Explanations of the data are given below.  Mean species
  lengths.  }
\usage{data(bimac)}
\format{
  A data frame with 45 observations on the following 8 variables.
  \describe{
    \item{node}{Labels for the nodes.}
    \item{species}{Species names for extant species.}
    \item{size}{Body size (head length in mm) of extant species.}
    \item{ancestor}{Ancestral node.}
    \item{time}{Time of node.}
    \item{OU.1}{a factor with levels \code{ns}}
    \item{OU.3}{a factor with levels \code{small} \code{medium} \code{large}}
    \item{OU.4}{a factor with levels \code{small} \code{medium} \code{large} \code{anc}}
    \item{OU.LP}{a factor with levels \code{small} \code{medium} \code{large}}
  }
}
\details{
  \item{Body size.}{We use the phenotypic data and phylogeny of Losos
  (1990), which employed the head lengths (of males) as a proxy for
  body size. In this group of lizards, head length correlates very
  strongly with snout-to-vent length and the cube root of mass, which
  are standard measures of body size. The data are head lengths in mm,
  note that we use the log of this value in analyses.}  

  \item{Tree topology}{The tree topology is encoded via two vectors:
  \code{ancestor} and \code{time}.  Each node of the phylogenetic tree
  has a corresponding row in the data frame, numbered from 1 to 45.
  The columns \code{ancestor} and \code{time} specify the phylogeny.
  The \code{ancestor} variable specifies the topology: it is a list
  indicating the ancestor of each node.  The root node has ancestor 0.
  The variable \code{time} specifies the temporal location of each node,
  with the root node being at time 0.}

  \item{Specifications of selective regimes.}{(Columns \code{OU.1},
    \code{OU.3}, \code{OU.4}, \code{OU.LP}).  These columns are
    factors, the levels of which correspond to the ``paintings'' of
    the respective adaptive regime hypotheses onto the phylogeny. Each
    selective regime is named (small, medium, large, etc.).  Put the
    corresponding name on each branch segment to indicate which
    selective regime it belongs to. Each column corresponds to a
    different painting of the selective regimes, and thus to a
    different hypothesis.  In this example, there are 3 alternative
    models (see Butler & King 2004): \code{OU.4} is 4-regime model,
    \code{OU.3} is 3-regime model (all ancestors are medium),
    \code{OU.LP} is linear parsimony model.}  }
\source{
  Butler, M.A. and A.A. King (2004) Phylogenetic comparative analysis: a
  modeling approach for adaptive evolution. American Naturalist, in press.
}
\references{
  Lazell, J. D. 1972. The anoles (Sauria: Iguanidae) of the Lesser
  Antilles. Bull. Mus. Comp. Zool., 143:1-115.

  Losos, J. B. 1990. A phylogenetic analysis of character displacement in
  Caribbean Anolis lizards. Evolution, 44:558-569.
}
\examples{
data(bimac)
attach(bimac)
tree.plot(node,ancestor,time,names=species,regimes=OU.LP)
brown.fit(log(size),node,ancestor,time/max(time))
hansen.fit(log(size),node,ancestor,time/max(time))
## a simple bootstrap:
sapply(hansen.dev(n=5,node,ancestor,time/max(time),alpha=0.19,sigma=0.22,theta=3.0),
       function(x)unlist(hansen.fit(x,node,ancestor,time/max(time))))
sapply(brown.dev(n=5,node,ancestor,time/max(time),sigma=0.22,theta=3),
       function(x)unlist(hansen.fit(x,node,ancestor,time/max(time))))
hansen.fit(log(size),node,ancestor,time/max(time),OU.3)
hansen.fit(log(size),node,ancestor,time/max(time),OU.4)
hansen.fit(log(size),node,ancestor,time/max(time),OU.LP)
detach(bimac)}
\author{Marguerite A. Butler <mabutler at utk dot edu> and Aaron A. King <kingaa at umich dot edu>}
\keyword{models}
