\name{isAvailableFuncCLP}
\alias{isAvailableFuncCLP}

\title{
  Check if some newer functions relying on newer CLP versions are available
}

\description{
  Checks if functions \code{setRowNameCLP}, \code{setColNameCLP}, 
  \code{writeMPSCLP}, and \code{modifyCoefficientCLP} are available.
}

\usage{
  isAvailableFuncCLP(funcname)
}

\arguments{
  \item{funcname}{
    The name of the function.
  }
}

\details{
  Some functions of the API require that the package was build against
  newer CLP versions with partly optional builds, otherwise they will be 
  dummy functions. These functions are \code{setRowNameCLP}, \code{setColNameCLP}, 
  \code{writeMpsCLP}, and \code{modifyCoefficientCLP}. This function 
  tests, if the respective function is available.
}

\value{
  TRUE, if the respective function is available, 
  FALSE, if it is not. If the funcname
  does not match one of the mentioned functions, it
  returns FALSE.
}

\author{
  Mayo Roettger <mayo.roettger@hhu.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

