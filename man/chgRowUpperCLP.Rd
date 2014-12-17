\name{chgRowUpperCLP}
\alias{chgRowUpperCLP}
\alias{Clp_chgRowUpper}

\title{
  Set/Change Row Upper Bounds
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_chgRowUpper}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
chgRowUpperCLP(lp, rub)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{rub}{
    Numeric vector containing the upper bounds of the
    rows of the model.
  }
}

\details{
  Interface to the C function \code{chgRowUpper} which calls
  the COIN-OR Clp function \code{Clp_chgRowUpper}.
}

\value{
  NULL
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

