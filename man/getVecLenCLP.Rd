\name{getVecLenCLP}
\alias{getVecLenCLP}
\alias{Clp_getVectorLengths}

\title{
  Retrieve the Number of Non Zero Elements per Column
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_getVectorLengths}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getVecLenCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getVecLen} which calls
  the COIN-OR Clp functions \code{Clp_numberColumns} and
  \code{Clp_getVectorLengths}.
}

\value{
  An integer vector containing the number of non zero
  elements per column.
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

