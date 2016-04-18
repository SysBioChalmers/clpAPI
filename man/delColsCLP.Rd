\name{delColsCLP}
\alias{delColsCLP}
\alias{Clp_deleteColumns}

\title{
  Delete Columns in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_deleteColumns}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  delColsCLP(lp, num, j)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{num}{
    Number of columns to delete.
  }
  \item{j}{
    Integer vector, containing the indices of columns to delete
    (the first column has index 0).
  }
}

\details{
  Interface to the C function \code{delCols} which calls
  the COIN-OR Clp function \code{Clp_deleteColumns}.
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

