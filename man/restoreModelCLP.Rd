\name{restoreModelCLP}
\alias{restoreModelCLP}
\alias{Clp_restoreModel}

\title{
  Restore model from file
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_restoreModel}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  restoreModelCLP(lp, fname)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{fname}{
    A filename.
  }
}

\details{
  Interface to the C function \code{restoreModel} which calls
  the COIN-OR Clp function \code{Clp_restoreModel}.
}

\value{
  Returns zero on success, otherwise non zero.
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

