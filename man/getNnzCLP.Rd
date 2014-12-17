\name{getNnzCLP}
\alias{getNnzCLP}
\alias{Clp_getElements}

\title{
  Retrieve the Non Zero Elements of the Constraint Matrix
  in Column Major Order.
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_getElements}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getNnzCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getNnz} which calls
  the COIN-OR Clp functions \code{Clp_getNumElements} and
  \code{Clp_getElements}.
}

\value{
  A numeric vector containing the non zero elements of the
  constraint matrix in column major order.
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

