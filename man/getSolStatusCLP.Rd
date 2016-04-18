\name{getSolStatusCLP}
\alias{getSolStatusCLP}
\alias{Clp_status}

\title{
  Retrieve the Solution Status
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_status}. Consult the COIN-OR Clp
  documentation for more detailed information.
}

\usage{
  getSolStatusCLP(lp)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
}

\details{
  Interface to the C function \code{getSolStatus} which calls
  the COIN-OR Clp function \code{Clp_status}.
}

\value{
  The solution status: 0: optimal, 1: primal infeasible,
  2: dual infeasible, 3: stopped on iterations etc,
  4: stopped due to errors.
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

