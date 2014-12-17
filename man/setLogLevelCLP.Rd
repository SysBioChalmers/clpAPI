\name{setLogLevelCLP}
\alias{setLogLevelCLP}
\alias{Clp_setLogLevel}

\title{
  Set the Amount of Output to STDOUT
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_setLogLevel}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  setLogLevelCLP(lp, amount)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{amount}{
    Log level flag: 0: nothing, 1: just final, 2: just
    factorizations, 3: as 2 plus a bit more, 4: verbose.
  }
}

\details{
  Interface to the C function \code{setLogLevel} which calls
  the COIN-OR Clp function \code{Clp_setLogLevel}.
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

