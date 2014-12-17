\name{dropNamesCLP}
\alias{dropNamesCLP}
\alias{Clp_dropNames}

\title{
  Drop Names in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_dropNames}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  dropNamesCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{dropNames} which calls
  the COIN-OR Clp function \code{Clp_dropNames}.
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

