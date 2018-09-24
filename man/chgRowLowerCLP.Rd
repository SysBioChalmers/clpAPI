\name{chgRowLowerCLP}
\alias{chgRowLowerCLP}
\alias{Clp_chgRowLower}

\title{
  Set/Change Row Lower Bounds
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_chgRowLower}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  chgRowLowerCLP(lp, rlb)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{rlb}{
    Numeric vector containing the lower bounds of the
    rows of the model.
  }
}

\details{
  Interface to the C function \code{chgColLower} which calls
  the COIN-OR Clp function \code{Clp_chgRowLower}.
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

