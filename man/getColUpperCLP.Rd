\name{getColUpperCLP}
\alias{getColUpperCLP}
\alias{Clp_columnUpper}

\title{
  Retrieve Column Upper Bounds
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_columnUpper}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  getColUpperCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getColLower} which calls
  the COIN-OR Clp functions \code{Clp_numberColumns} and
  \code{Clp_columnUpper}.
}

\value{
  The upper bounds of the models columns (the corresponding
  structural variables) are returned.
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

