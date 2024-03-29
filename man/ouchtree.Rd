% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ouchtree.R, R/ape2ouch.R
\name{ouchtree}
\alias{ouchtree}
\alias{ouchtree-class}
\alias{ape2ouch}
\title{Phylogenetic tree object in \pkg{ouch} format}
\usage{
ouchtree(nodes, ancestors, times, labels = as.character(nodes))

ape2ouch(tree, scale = TRUE, branch.lengths = tree$edge.length)
}
\arguments{
\item{nodes}{A character vector giving the name of each node.
These are used internally and must be unique.}

\item{ancestors}{Specification of the topology of the phylogenetic tree.
This is in the form of a character vector specifying the name
(as given in the \code{nodes} argument)
of the immediate ancestor of each node.
In particular, the i-th name is that of the ancestor of the i-th node.
The root node is distinguished by having no ancestor (i.e., \code{NA}).}

\item{times}{A vector of nonnegative numbers, one per node in the tree,
specifying the time at which each node is located.
Time should be increasing from the root node to the terminal twigs.}

\item{labels}{Optional vector of node labels.
These will be used in plots to label nodes.
It is not necessary that these be unique.}

\item{tree}{a tree of class \link[ape:read.tree]{ape::phylo}.}

\item{scale}{optional.
If \code{scale=TRUE}, the tree's depth will be scaled to 1.
If \code{scale} is a number, then the branch lengths will be scaled by this number.}

\item{branch.lengths}{optional vector of branch lengths.}
}
\description{
\code{ouchtree} constructs a representation of a phylogenetic tree.

\code{ape2ouch} translates \pkg{ape}'s \code{phylo} representation of a phylogenetic tree into \pkg{ouch}'s \code{ouchtree} representation.
Optionally, the user can adjust the branch lengths while preserving the topology.
}
\details{
\code{ouchtree()} creates an \code{ouchtree} object given information on the phylogeny's topology and node times.
An \code{ouchtree} object also (optionally) holds names of taxa for display purposes.
}
\examples{
tree <- with(
  bimac,
  ouchtree(nodes=node,ancestors=ancestor,times=time,labels=spcode)
)
tree

plot(tree)
plot(tree, node.names=TRUE)    # display node names

}
\seealso{
Other phylogenetic comparative models: 
\code{\link{brown}()},
\code{\link{hansen}()},
\code{\link{ouch-package}},
\code{\link{paint}()}
}
\author{
Aaron A. King

A. A. King, D. Ackerly
}
\concept{phylogenetic comparative models}
\keyword{models}
