\name{primalCLP}
\alias{primalCLP}
\alias{Clp_primal}

\title{
  Solve LP Problem with the Primal Simplex Method
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_primal}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  primalCLP(lp, ifValP = 0)
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
  Interface to the C function \code{primal} which calls
  the COIN-OR Clp function \code{Clp_primal}.
}

\value{
  A return code.
}

\references{
   The COIN-OR Clp home page at
   \url{http://www.coin-or.org/projects/Clp.xml}
}

\author{
  Gabriel Gelius-Dietrich

  Maintainer: Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>
}

\keyword{ optimize }

