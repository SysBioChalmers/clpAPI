\name{printModelCLP}
\alias{printModelCLP}
\alias{Clp_printModel}

\title{
  Print the Model to STDOUT
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_printModel}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  printModelCLP(lp, prefix = "CLPmodel")
}


\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{prefix}{
    A character string containing a name for the model.
  }
}

\details{
  Interface to the C function \code{printModel} which calls
  the COIN-OR Clp function \code{Clp_printModel}.
}

\value{
  NULL
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

