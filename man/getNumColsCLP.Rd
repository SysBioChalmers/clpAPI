\name{getNumColsCLP}
\alias{getNumColsCLP}
\alias{Clp_numberColumns}

\title{
  Retrieve the Current Number of Columns in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_numberColumns}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getNumColsCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getNumCols} which calls
  the COIN-OR Clp function \code{Clp_numberColumns}.
}

\value{
  The current number of columns in the model.
}

\references{
   The COIN-OR Clp home page at
   \url{https://github.com/coin-or/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

