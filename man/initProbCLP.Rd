\name{initProbCLP}
\alias{initProbCLP}
\alias{Clp_newModel}

\title{
  Create a COIN-OR Clp Problem Object
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_newModel}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  initProbCLP(ptrtype = "clp_prob")
}

\arguments{
  \item{ptrtype}{
    A name for the pointer to a COIN-OR Clp problem object.
  }
}

\details{
  Interface to the C function \code{initProb} which calls
  the COIN-OR Clp function \code{Clp_newModel}.
}

\value{
  An instance of class \code{"\linkS4class{clpPtr}"}.
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

