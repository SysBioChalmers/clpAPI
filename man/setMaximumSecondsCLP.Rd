\name{setMaximumSecondsCLP}
\alias{setMaximumSecondsCLP}
\alias{Clp_setMaximumSeconds}

\title{
  Set the Maximum Time in Seconds
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_setMaximumSeconds}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  setMaximumSecondsCLP(lp, seconds)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{seconds}{
    Maximum duration in seconds
  }
}

\details{
  Interface to the C function \code{setMaximumSeconds} which calls
  the COIN-OR Clp function \code{Clp_setMaximumSeconds}.
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

