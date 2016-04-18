\name{probNameCLP}
\alias{probNameCLP}
\alias{Clp_problemName}

\title{
  Set Problem Name
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_problemName}. Consult the COIN-OR Clp documentation
  for more detailed information.

}

\usage{
  probNameCLP(lp, pname)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{pname}{
    A single character string containing the problem name.
  }
}

\details{
  Interface to the C function \code{probName} which calls
  the COIN-OR Clp function \code{Clp_problemName}.
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
