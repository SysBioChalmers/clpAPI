\name{getObjValCLP}
\alias{getObjValCLP}
\alias{Clp_objectiveValue}

\title{
  Retrieve the Value of the Objective Function After Optimization
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_objectiveValue}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getObjValCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getObjVal} which calls
  the COIN-OR Clp function \code{Clp_objectiveValue}.
}

\value{
  Returns the value of the objective function after optimization.
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

