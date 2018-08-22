\name{getRowUpperCLP}
\alias{getRowUpperCLP}
\alias{Clp_rowUpper}

\title{
  Retrieve Row Upper Bound
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_rowUpper}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  getRowUpperCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getRowUpper} which calls
  the COIN-OR Clp functions \code{Clp_numberRows} and
  \code{Clp_rowUpper}.
}

\value{
  The upper bounds of the models rows are returned.
}

\references{
   The COIN-OR Clp home page at
   \url{https://projects.coin-or.org/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

