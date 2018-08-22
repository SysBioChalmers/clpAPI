\name{solveInitialBarrierCLP}
\alias{solveInitialBarrierCLP}
\alias{Clp_initialBarrierSolve}

\title{
  Solve LP Problem with the Initial Barrier Method
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_initialBarrierSolve}. Consult the   COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  solveInitialBarrierCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{solveInitialBarrier}
  which calls the COIN-OR Clp function
  \code{Clp_initialBarrierSolve}.
}

\value{
  A return code.
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

