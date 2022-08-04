\name{setRowNameCLP}
\alias{setRowNameCLP}
\alias{Clp_setRowName}

\title{
  Set row name
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_setRowName}. Consult the COIN-OR Clp documentation
  for more detailed information. This function is only available 
  in newer COIN-OR Clp versions and/or specific builds. See 
  \code{isAvailableFuncCLP()} to check for availability.
}

\usage{
  setRowNameCLP(lp, i, rname)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }

  \item{i}{
    Row index.
  }
  
  \item{rname}{
    A single character string containing the row name.
  }
}

\details{
  Interface to the C function \code{setRowName} which calls
  the COIN-OR Clp function \code{Clp_setRowName}.
}

\value{
  NULL
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

