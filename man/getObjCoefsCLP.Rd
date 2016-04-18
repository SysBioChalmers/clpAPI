\name{getObjCoefsCLP}
\alias{getObjCoefsCLP}
\alias{Clp_objective}

\title{
  Retrieve Objective Coefficients
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_objective}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  getObjCoefsCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getObjCoefs} which calls
  the COIN-OR Clp functions \code{Clp_numberColumns} and
  \code{Clp_objective}.
}

\value{
  A numeric vector containing the objective coefficients.
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

