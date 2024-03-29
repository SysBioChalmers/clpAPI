\name{getScaleFlagCLP}
\alias{getScaleFlagCLP}
\alias{Clp_scalingFlag}

\title{
  Retrieve the Scale Flag
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_scalingFlag}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getScaleFlagCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getScaleFlag} which calls
  the COIN-OR Clp function \code{Clp_scalingFlag}.
}

\value{
  Returns the scaling flag: 0: off, 1: equilibrium, 2: geometric,
  3: auto, 4: dynamic (later - maybe not implemented in CLP?).
}

\references{
   The COIN-OR Clp home page at
   \url{https://github.com/coin-or/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

