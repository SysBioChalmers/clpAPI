\name{addRowsCLP}
\alias{addRowsCLP}
\alias{Clp_addRows}

\title{
  Add Rows
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_addRows}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  addRowsCLP(lp, nrows, lb, ub, rowst, cols, val)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{nrows}{
    Number of rows to add.
  }
  \item{lb}{
    Lower bounds of the new rows.
  }
  \item{ub}{
    Upper bounds of the new rows.
  }
  \item{rowst}{
    Vector containing the starting indices of new rows (Arguments \code{cols}
    and \code{val} must be in row major order). The first element of
    \code{rowst} must be \code{0}, the last element must be
    \code{length(val)+1}.
  }
  \item{cols}{
    Column indices of the non zero elements in the new rows.
  }
  \item{val}{
    Numerical values of the new non zero elements.
  }
}

\details{
  Interface to the C function \code{addRows} which calls
  the COIN-OR Clp function \code{Clp_addRows}.
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
