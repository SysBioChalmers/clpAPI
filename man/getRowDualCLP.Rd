\name{getRowDualCLP}
\alias{getRowDualCLP}
\alias{Clp_dualRowSolution}

\title{
  Retrieve all Row Dual Values
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_dualRowSolution}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getRowDualCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getRowDual} which calls
  the COIN-OR Clp functions \code{Clp_numberRows} and
  \code{Clp_dualRowSolution}.
}

\value{
  Returns all dual values of the auxiliary variables as a numeric vector.
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

