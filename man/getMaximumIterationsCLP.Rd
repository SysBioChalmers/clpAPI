\name{getMaximumIterationsCLP}
\alias{getMaximumIterationsCLP}
\alias{maximumIterations}

\title{
  Returns Maximum Number of Iterations
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{maximumIterations}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getMaximumIterationsCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getMaximumIterations} which calls
  the COIN-OR Clp function \code{maximumIterations}.
}

\value{
  Maximum number of iterations
}

\references{
   The COIN-OR Clp home page at
   \url{https://projects.coin-or.org/Clp}
}

\author{
  C. Jonathan Fritzemeier <clausjonathan.fritzemeier@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

