\name{modifyCoefficientCLP}
\alias{modifyCoefficientCLP}
\alias{Clp_modifyCoefficient}

\title{
  Change matrix coefficients
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_modifyCoefficient}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  modifyCoefficientCLP(lp, i, j, el, keepZero = TRUE)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{i}{
    Row index.
  }
  \item{j}{
    Column index.
  }
  \item{el}{
    Coefficient to set.
  }
  \item{keepZero}{
    If set to \code{TRUE}, keep zeroes.
  }
}

\details{
  Interface to the C function \code{modifyCoefficient} which calls
  the COIN-OR Clp function \code{Clp_modifyCoefficient}.
}

\value{
  NULL
}

\references{
   The COIN-OR Clp home page at
   \url{https://projects.coin-or.org/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

