\name{chgObjCoefsCLP}
\alias{chgObjCoefsCLP}
\alias{Clp_chgObjCoefficients}

\title{
  Set/Change Objective Coefficients
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_chgObjCoefficients}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  chgObjCoefsCLP(lp, objCoef)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{objCoef}{
    Numeric vector containing the objective coefficients
    of the model.
  }
}

\details{
  Interface to the C function \code{chgObjCoefs} which calls
  the COIN-OR Clp function \code{Clp_chgObjCoefficients}.
}

\value{
  NULL
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

