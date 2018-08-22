\name{getRowPrimCLP}
\alias{getRowPrimCLP}
\alias{Clp_primalRowSolution}

\title{
  Retrieve all Row Primal Values
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_primalRowSolution}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getRowPrimCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getRowPrim} which calls
  the COIN-OR Clp functions \code{Clp_numberRows} and
  \code{Clp_primalRowSolution}.
}

\value{
  Returns all primal values of the auxiliary variables as a numeric vector.
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

