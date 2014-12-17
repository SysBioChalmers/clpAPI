\name{solveInitialBarrierNoCrossCLP}
\alias{solveInitialBarrierNoCrossCLP}
\alias{Clp_initialBarrierNoCrossSolve}

\title{
  Solve LP Problem with the Initial Barrier Method
  (no Crossover)
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_initialBarrierNoCrossSolve}. Consult the
  COIN-OR Clp documentation for more detailed information.
}

\usage{
  solveInitialBarrierNoCrossCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{solveInitialBarrierNoCross}
  which calls the COIN-OR Clp function
  \code{Clp_initialBarrierNoCrossSolve}.
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

