\name{getIndCLP}
\alias{getIndCLP}
\alias{Clp_getIndices}

\title{
  Retrieve Row Indices of the Non Zero Elements in the Constraint Matrix
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_getIndices}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getIndCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getInd} which calls
  the COIN-OR Clp functions \code{Clp_numberColumns} and
  \code{Clp_getIndices}.
}

\value{
  An integer vector containing the row Indices of the
  non zero elements in the constraint matrix.
}

\references{
   The COIN-OR Clp home page at
   \url{https://github.com/coin-or/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

