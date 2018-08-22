\name{getColLowerCLP}
\alias{getColLowerCLP}
\alias{Clp_columnLower}

\title{
  Retrieve Column Lower Bound
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_columnLower}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  getColLowerCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getColLower} which calls
  the COIN-OR Clp functions \code{Clp_numberColumns} and
  \code{Clp_columnLower}.
}

\value{
  The lower bounds of the models columns (the corresponding
  structural variables) are returned.
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

