\name{delRowsCLP}
\alias{delRowsCLP}
\alias{Clp_deleteRows}

\title{
  Delete Rows in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_deleteRows}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  delRowsCLP(lp, num, i)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{num}{
    Number of rows to delete.
  }
  \item{i}{
    Integer vector, containing the indices of rows to delete
    (the first row has index 0).
  }
}

\details{
  Interface to the C function \code{delRows} which calls
  the COIN-OR Clp function \code{Clp_deleteRows}.
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

