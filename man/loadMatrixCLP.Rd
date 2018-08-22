\name{loadMatrixCLP}
\alias{loadMatrixCLP}

\title{
  Load Constraint Matrix
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_loadProblem}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  loadMatrixCLP(lp, ncols, nrows, ia, ja, ra)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{ncols}{
    Number of Columns.
  }
  \item{nrows}{
    Number of Rows.
  }
  \item{ia}{
    Row indices in the constraint matrix.
  }
  \item{ja}{
    Column starts in constraint matrix.
  }
  \item{ra}{
    Non zero elements of the constraint matrix.
  }
}

\details{
  Interface to the C function \code{loadMatrix} which calls
  the COIN-OR Clp function \code{Clp_loadProblem}.
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

