\name{chgColUpperCLP}
\alias{chgColUpperCLP}
\alias{Clp_chgColumnUpper}

\title{
  Set/Change Column Upper Bounds
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_chgColumnUpper}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  chgColUpperCLP(lp, ub)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{ub}{
    Numeric vector containing the upper bounds of the
    columns of the model.
  }
}

\details{
  Interface to the C function \code{chgColUpper} which calls
  the COIN-OR Clp function \code{Clp_chgColumnUpper}.
}

\value{
  NULL
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

