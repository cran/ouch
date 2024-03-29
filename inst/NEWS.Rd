\name{NEWS}
\title{News for Package 'ouch'}
\section{Changes in ouch version 2.19}{
  \itemize{
    \item There is a new \code{as.data.frame} method.
  }
}
\section{Changes in ouch version 2.17}{
  \itemize{
    \item Both left and right plot margins can be adjusted to allow, for example, for long taxon names or long legend labels.
    \item There is now an example showing how to plot a tree with long taxon names.
    \item Error messages and examples have been improved somewhat.
    \item The manual and help pages have been streamlined.
  }
}
\section{Changes in ouch version 2.16}{
  \itemize{
    \item The arguments to the various \code{plot} methods have changed.
    There is now much finer control over the plotting.
    In particular, one can specify a different palette of colors, adjust the width of the margin in which taxon labels are printed, and pass options separately to the base-graphics \code{plot}, \code{text}, and \code{legend} commands under0lying the \pkg{ouch} \code{plot} methods.
    \item A layout bug in \code{plot} has been fixed.
  }
}
\section{Changes in ouch version 2.15}{
  \itemize{
    \item The \code{plot} function now spaces the terminal taxa evenly.
    \item The documentation has received an overhaul.
  }
}
\section{Changes in ouch version 2.14}{
  \itemize{
    \item Depends on R>=3.6.
    This change necessitated by modifications to the default RNG with R version 3.6.
  }
}
\section{Changes in ouch version 2.13}{
  \itemize{
    \item Changes to the organization of the package.
  }
}
\section{Changes in ouch version 2.12}{
  \itemize{
    \item Upgrade the package license to GPL version 3.
  }
}
\section{Changes in ouch version 2.11}{
  \itemize{
    \item Improvements to the underlying C code, which now registers callable functions.
  }
}
\section{Changes in ouch version 2.10}{
  \itemize{
    \item In \code{glssoln}, when the Choleski decomposition fails, a warning is issued and \code{NA}s are returned.
    In previous versions, an error was thrown.
  }
}
\section{Changes in ouch version 2.8}{
  \itemize{
    \item 'as(tree,"data.frame")' now returns a data-frame with rownames identical to the node-names (tree@nodes).
    \item Functions that have been deprecated since 2.1-1 have been removed.
      The argument 'alpha' that was deprecated in favor of 'sqrt.alpha' in version 2.6-1 has been removed.
  }
}
\section{Changes in ouch version 2.7}{
  \itemize{
    \item Fixed typo in the documentation of 'hansen'.  Thanks to Roland Sookias for bringing this to my attention.
  }
}
\section{Changes in ouch version 2.6}{
  \itemize{
    \item In the 'hansen' function, the use of the argument 'alpha' was confusing.  
      	This vector parameterizes the "selection" portion of the OU model, which since version 2 has been assumed to be a positive-definite matrix (corresponding to stabilizing selection).
      	In particular, it contains the entries of the lower-triangular matrix square root of this positive definite matrix.
      To make the name more mnemonic, this argument has been renamed 'sqrt.alpha'.
      This is a superficial change meant to improve intelligibility: the algorithm itself has not changed.
      To prevent loss of backward compatibility, use of argument 'alpha' will continue to work, though it will generate a warning.
      In some future release, use of 'alpha' as an argument to 'hansen' will no longer generate a warning but will instead be passed through to the control list for 'optim' if method='optim' and will otherwise generate an error.
  }
}
\section{Changes in ouch version 2.1}{
  \itemize{
    \item Version 2 of ouch represents a major revision.  
      	ouch's capabilities have been substantially extended.
    \item support is now provided for multivariate characters.  
          Inference based on multivariate Ornstein-Uhlenbeck and multivariate Brownian motion models is now provided.
    \item nonultrametric trees are supported.
    \item the root character state need no longer be estimated.  
          Instead, the root state is assumed to be distributed according to the stationary distribution of the OU process.
    \item convenience functions for parametric bootstrap confidence intervals are provided.
    \item the new function 'paint' simplifies the specification of adaptive regimes in complex phylogenies.
    \item the new function 'ape2ouch' translates trees in the 'phylo' format of the ape package into 'ouchtree' format (thanks to Dan Ackerly).
    \item The new package is a complete rewrite, using S4 classes.  
      	The old functionality from ouch v. 1 remains available, but is now deprecated.
  }
}
