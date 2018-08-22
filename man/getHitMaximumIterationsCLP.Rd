\name{getHitMaximumIterationsCLP}
\alias{getHitMaximumIterationsCLP}
\alias{Clp_hitMaximumIterations}

\title{
  Returns True if Hit Maximum Iterations (or Time)
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_hitMaximumIterations}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getHitMaximumIterationsCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getHitMaximumIterations} which calls
  the COIN-OR Clp function \code{Clp_hitMaximumIterations}.
}

\value{
  True if hit maximum iterations (or time) 
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

