\name{addColsCLP}
\alias{addColsCLP}
\alias{Clp_addColumns}

\title{
  Add Columns
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_addColumns}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  addColsCLP(lp, ncols, lb, ub, obj, colst, rows, val)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{ncols}{
    Number of columns to add.
  }
  \item{lb}{
    Lower bounds of the new columns.
  }
  \item{ub}{
    Upper bounds of the new columns.
  }
  \item{obj}{
    Objective coeficients of the new columns.
  }
  \item{colst}{
    Vector containing the starting indices of new rows (Arguments \code{rows}
    and \code{val} must be in column major order). The first element of
    \code{colst} must be \code{0}, the last element must be
    \code{length(val)+1}.
  }
  \item{rows}{
    Row indices of the non zero elements in the new columns.
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
