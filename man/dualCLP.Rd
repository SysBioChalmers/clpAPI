\name{dualCLP}
\alias{dualCLP}
\alias{Clp_dual}

\title{
  Solve LP Problem with the Dual Simplex Method
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_dual}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  dualCLP(lp, ifValP = 0)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{ifValP}{
    An integer value.
  }
}

\details{
  Interface to the C function \code{dual} which calls
  the COIN-OR Clp function \code{Clp_dual}.
}

\value{
  A return code.
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

