\name{resizeCLP}
\alias{resizeCLP}
\alias{Clp_resize}

\title{
  Resize a Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_resize}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  resizeCLP(lp, nrows, ncols)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{nrows}{
    Number of rows.
  }
  \item{ncols}{
    Number of columns.
  }
}

\details{
  Interface to the C function \code{resize} which calls
  the COIN-OR Clp function \code{Clp_resize}.
  
  The function \code{resize} can produce a larger model.
  If the current number of rows and columns is \emph{n}
  and \emph{m} respectively and you set nrows to \emph{i}
  and ncols to \emph{j}, the new number of rows and columns
  will be \emph{i} and \emph{j}. It is not possible to
  scale down the model. In order to delete rows or columns,
  use \code{\link{delRowsCLP}} or \code{\link{delColsCLP}}.
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


\seealso{
  \code{\link{delRowsCLP}} and \code{\link{delColsCLP}}.
}

\keyword{ optimize }

