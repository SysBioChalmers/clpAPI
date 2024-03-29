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
   \url{https://github.com/coin-or/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

