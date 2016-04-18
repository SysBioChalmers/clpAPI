\name{getNumNnzCLP}
\alias{getNumNnzCLP}
\alias{Clp_getNumElements}

\title{
  Retrieve the Current Number of Non Zero Elements in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_getNumElements}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getNumNnzCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getNumNnz} which calls
  the COIN-OR Clp function \code{Clp_getNumElements}.
}

\value{
  Returns the current number of non zero elements in the model.
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

