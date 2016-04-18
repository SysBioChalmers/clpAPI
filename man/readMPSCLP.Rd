\name{readMPSCLP}
\alias{readMPSCLP}
\alias{Clp_readMps}

\title{
  Read Problem in (Free) MPS Format
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_readMps}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  readMPSCLP(lp, fname, keepNames = TRUE, ignoreErrors = FALSE)
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
  \item{keepNames}{
    Boolean, keep variable names.
  }
  \item{ignoreErrors}{
    If set to \code{TRUE}, errors will be ignored.
  }
}

\details{
  Interface to the C function \code{readMPS} which calls
  the COIN-OR Clp function \code{Clp_readMps}.
}

\value{
  Returns zero on success, otherwise non zero.
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

