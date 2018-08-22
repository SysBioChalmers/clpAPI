\name{solveInitialCLP}
\alias{solveInitialCLP}
\alias{Clp_initialSolve}

\title{
  Solve LP Problem with a General Solve Algorithm
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_initialSolve}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  solveInitialCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{solveInitial}
  which calls the COIN-OR Clp function
  \code{Clp_initialSolve}.
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

