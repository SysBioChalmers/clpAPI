\name{getLogLevelCLP}
\alias{getLogLevelCLP}
\alias{Clp_logLevel}

\title{
  Retrieve the Log Level Flag
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_logLevel}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getLogLevelCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getLogLevel} which calls
  the COIN-OR Clp function \code{Clp_logLevel}.
}

\value{
  Returns the log level flag: 0: nothing, 1: just final, 2: just
  factorizations, 3: as 2 plus a bit more, 4: verbose.
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

