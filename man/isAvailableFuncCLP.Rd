\name{isAvailableFuncCLP}
\alias{isAvailableFuncCLP}

\title{
  Check if some newer functions relying on CLP version >= 1.17.2 are available
}

\description{
  Checks if functions like \code{setRowNameCLP}, \code{setColNameCLP}, 
  \code{writeMpsCLP}, and \code{modifyCoefficientCLP} are available.
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
  CLP version >= 1.17.2, otherwise they will be dummy functions. 
  These functions are \code{setRowNameCLP}, \code{setColNameCLP}, 
  \code{writeMpsCLP}, and \code{modifyCoefficientCLP}. This function 
  tests, if these functions are available.
}

\value{
  TRUE, if the respective function is available, 
  returns FALSE, if it is not. If the funcname
  does not match one of the mentioned functions, it
  returns FALSE.
}

\author{
  Mayo Roettger <mayo.roettger@hhu.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\keyword{ optimize }

