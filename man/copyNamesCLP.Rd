\name{copyNamesCLP}
\alias{copyNamesCLP}
\alias{Clp_copyNames}

\title{
  Copy Column and Row Names in the Model
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_copyNames}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  copyNamesCLP(lp, cnames, rnames)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{cnames}{
    Character vector, containing the column names, must not be longer than the
    number of columns in the model.
  }
  \item{rnames}{
    Character vector, containing the row names, must not be longer than the
    number of rows in the model.
  }
}

\details{
  Interface to the C function \code{copyNames} which calls
  the COIN-OR Clp function \code{Clp_copyNames}.
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

