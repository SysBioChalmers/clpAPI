\name{loadProblemCLP}
\alias{loadProblemCLP}
\alias{Clp_loadProblem}

\title{
  Load Problem Data
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_loadProblem}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  loadProblemCLP(lp, ncols, nrows, ia, ja, ra,
                 lb = NULL, ub = NULL, obj_coef = NULL,
                 rlb = NULL, rub = NULL)
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
  \item{lb}{
    Column lower bounds.
  }
  \item{ub}{
    Column upper bounds.
  }
  \item{obj_coef}{
    Objective coefficients.
  }
  \item{rlb}{
    Row lower bounds.
  }
  \item{rub}{
    Row upper bounds.
  }
}

\details{
  Interface to the C function \code{loadProblem} which calls
  the COIN-OR Clp function \code{Clp_loadProblem}.
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

