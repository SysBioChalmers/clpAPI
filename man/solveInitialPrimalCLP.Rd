\name{solveInitialPrimalCLP}
\alias{solveInitialPrimalCLP}
\alias{Clp_initialPrimalSolve}

\title{
  Solve LP Problem with the Initial Primal Simplex Method
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_initialPrimalSolve}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  solveInitialPrimalCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{solveInitialPrimal}
  which calls the COIN-OR Clp function
  \code{Clp_initialPrimalSolve}.
}

\value{
  A return code.
}

\references{
   The COIN-OR Clp home page at
   \url{http://www.coin-or.org/projects/Clp.xml}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Claus Jonathan Fritzemeier <clausjonathan.fritzemeier@uni-duesseldorf.de>
}


\keyword{ optimize }

