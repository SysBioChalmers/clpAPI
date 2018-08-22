\name{setMaximumIterationsCLP}
\alias{setMaximumIterationsCLP}
\alias{Clp_setMaximumIterations}

\title{
  Set the Maximum Number of Iterations
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_setMaximumIterations}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  setMaximumIterationsCLP(lp, iterations)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{iterations}{
    Number of iterations
  }
}

\details{
  Interface to the C function \code{setMaximumIterations} which calls
  the COIN-OR Clp function \code{Clp_setMaximumIterations}.
}

\value{
  NULL
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

