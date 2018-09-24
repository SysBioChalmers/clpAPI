\name{getMaximumSecondsCLP}
\alias{getMaximumSecondsCLP}
\alias{Clp_maximumSeconds}

\title{
  Maximum Time in Seconds (from when Set called)
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_maximumSeconds}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getMaximumSecondsCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getMaximumSeconds} which calls
  the COIN-OR Clp function \code{Clp_maximumSeconds}.
}

\value{
  Maximum time in seconds (from when set called)
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

