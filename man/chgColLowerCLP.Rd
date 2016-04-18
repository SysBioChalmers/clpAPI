\name{chgColLowerCLP}
\alias{chgColLowerCLP}
\alias{Clp_chgColumnLower}

\title{
  Set/Change Column Lower Bounds
}

\description{
  Low level interface function to the COIN-OR Clp function
  \code{Clp_chgColumnLower}. Consult the COIN-OR Clp documentation
  for more detailed information.
}

\usage{
  chgColLowerCLP(lp, lb)
}

\arguments{
  \item{lp}{
    An object of class \code{"\linkS4class{clpPtr}"} as returned by
    \code{\link{initProbCLP}}. This is basically a pointer to a COIN-OR Clp
    problem object.
  }
  \item{lb}{
    Numeric vector containing the lower bounds of the
    columns of the model.
  }
}

\details{
  Interface to the C function \code{chgColLower} which calls
  the COIN-OR Clp function \code{Clp_chgColumnLower}.
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

