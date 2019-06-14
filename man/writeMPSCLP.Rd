\name{writeMPSCLP}
\alias{writeMPSCLP}
\alias{Clp_writeMps}

\title{
  Write an MPS Format file to the given filename
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_writeMps}. Consult the COIN-OR Clp documentation
  for more detailed information. This function is only available,
  if you are using COIN-OR Clp version >= 1.17.2.
}

\usage{
  writeMPSCLP(lp, fname, formatType = 0, numberAcross = 1, objSense = 1)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  
  \item{fname}{
    A filename.
  }
  
  \item{formatType}{
    Integer value: 0 = normal, 1 = extra or 2 = hex.
  }
  
  \item{numberAcross}{
    Number across is 1 or 2.
  }
  
  \item{objSense}{
    Use objSense = -1 to flip the objective funtion around.
  }
}

\details{
  Interface to the C function \code{writeMps} which calls
  the COIN-OR Clp function \code{Clp_writeMps}.
}

\value{
  Returns zero on success, otherwise non zero.
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

