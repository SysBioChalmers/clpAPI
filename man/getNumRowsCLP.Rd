\name{getNumRowsCLP}
\alias{getNumRowsCLP}
\alias{Clp_numberRows}

\title{
  Retrieve the Current Number of Rows in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_numberRows}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getNumRowsCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getNumRows} which calls
  the COIN-OR Clp function \code{Clp_numberRows}.
}

\value{
  The current number of rows in the model.
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

