\name{scaleModelCLP}
\alias{scaleModelCLP}
\alias{Clp_scaling}

\title{
  Set/Unset the Scaling Flag (Method)
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_scaling}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  scaleModelCLP(lp, mode)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{mode}{
    Scaling flag: 0: off, 1: equilibrium, 2: geometric,
    3: auto, 4: dynamic (later - maybe not implemented
    in CLP?).
  }
}

\details{
  Interface to the C function \code{scaleModel} which calls
  the COIN-OR Clp function \code{Clp_scaling}.
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

