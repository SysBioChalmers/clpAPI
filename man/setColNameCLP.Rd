\name{setColNameCLP}
\alias{setColNameCLP}
\alias{Clp_setColumnName}

\title{
  Set column name
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_setColumnName}. Consult the COIN-OR Clp documentation
  for more detailed information. This function is only available,
  if you are using COIN-OR Clp version >= 1.17.2.
}

\usage{
  setColNameCLP(lp, j, cname)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  
  \item{j}{
    Column index.
  }
  
  \item{cname}{
    A single character string containing the column name.
  }
}

\details{
  Interface to the C function \code{setColName} which calls
  the COIN-OR Clp function \code{Clp_setColumnName}.
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

