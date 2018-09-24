\name{lengthNamesCLP}
\alias{lengthNamesCLP}
\alias{Clp_lengthNames}

\title{
  Length of the Names in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_lengthNames}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  lengthNamesCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{lengthNames} which calls
  the COIN-OR Clp function \code{Clp_lengthNames}.
}

\value{
  Number of characters of the longest name in the Model.
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

