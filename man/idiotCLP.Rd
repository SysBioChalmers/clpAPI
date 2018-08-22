\name{idiotCLP}
\alias{idiotCLP}
\alias{Clp_idiot}

\title{
  Solve LP Problem with the idiot Code
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_idiot}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  idiotCLP(lp, thd = 0)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{thd}{
    An integer value.
  }
}

\details{
  Interface to the C function \code{idiot} which calls
  the COIN-OR Clp function \code{Clp_idiot}.
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

