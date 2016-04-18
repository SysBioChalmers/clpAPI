\name{clpAPI-package}
\alias{clpAPI-package}
\alias{clpAPI}

\docType{package}

\title{
  R Interface to C API of COIN-OR Clp
}

\description{
  A low level interface to COIN-OR Clp (COIN Linear Program code).
}

\details{
  The package \code{clpAPI} provides access to the callable library
  of COIN-OR Clp from within R.
}

\references{
   The COIN-OR Clp home page at
   \url{http://www.coin-or.org/projects/Clp.xml}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Claus Jonathan Fritzemeier <clausjonathan.fritzemeier@uni-duesseldorf.de>
}


\keyword{ package }
\keyword{ optimize}

\examples{
# load package
library(clpAPI)

# preparing the model
lp <- initProbCLP()

nrows  <- 5
ncols  <- 8

# objective function
obj    <- c(1, 0, 0, 0, 2, 0, 0, -1)

# upper and lower bounds of the rows
rlower <- c(2.5, -1000, 4, 1.8, 3)
rupper <- c(1000, 2.1, 4, 5, 15)

# upper and lower bounds of the columns
clower <- c(2.5, 0, 0, 0, 0.5, 0, 0, 0)
cupper <- c(1000, 4.1, 1, 1, 4, 1000, 1000, 4.3)

# constraint matrix
ia <- c(0, 4, 0, 1, 1, 2, 0, 3, 0, 4, 2, 3, 0, 4)
ja <- c(0, 2, 4, 6, 8, 10, 11, 12, 14)
ar <- c(3.0, 5.6, 1.0, 2.0, 1.1, 1.0, -2.0, 2.8,
        -1.0, 1.0, 1.0, -1.2, -1.0, 1.9)

# direction of optimization
setObjDirCLP(lp, 1)

# load problem data
loadProblemCLP(lp, ncols, nrows, ia, ja, ar,
               clower, cupper, obj, rlower, rupper)

# solve lp problem
solveInitialCLP(lp)

# retrieve the results
getSolStatusCLP(lp)
getObjValCLP(lp)
getColPrimCLP(lp)

# remove problem object
delProbCLP(lp)
}
