\name{getColPrimCLP}
\alias{getColPrimCLP}
\alias{Clp_primalColumnSolution}

\title{
  Retrieve all Column Primal Values
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_primalColumnSolution}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getColPrimCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getColPrim} which calls
  the COIN-OR Clp functions \code{Clp_numberColumns} and
  \code{Clp_primalColumnSolution}.
}

\value{
  Returns all primal values of the stuctural variables as a numeric vector.
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

