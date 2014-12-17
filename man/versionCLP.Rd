\name{versionCLP}
\alias{versionCLP}
\alias{CLP_VERSION}

\title{
  Determine COIN-OR Clp Callable Library Version
}

\description{
  Low level interface function to the COIN-OR Clp constant \code{CLP_VERSION}.
   Consult the COIN-OR Clp documentation for more detailed information.
}

\usage{
  versionCLP()
}

\details{
  Interface to the C function \code{version} which returns the COIN-OR Clp
  version number.
}

\value{
  Returns a single character value containing the COIN-OR Clp version number.
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
