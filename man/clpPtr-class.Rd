\name{clpPtr-class}

\Rdversion{1.1}
\docType{class}

\alias{clpPtr-class}
\alias{clpPtr}
\alias{isCLPpointer,clpPtr-method}
\alias{isNULLpointerCLP,clpPtr-method}
\alias{clpPointer,clpPtr-method}
\alias{clpPtrType,clpPtr-method}
\alias{clpPtrType<-,clpPtr-method}
\alias{isCLPpointer}
\alias{isNULLpointerCLP}
\alias{clpPointer}
\alias{clpPtrType}
\alias{clpPtrType<-}

\title{Class \code{"clpPtr"}}

\description{
  Structure of the class \code{"clpPtr"}. Objects of that class are
  used to hold pointers to C structures used by COIN-OR Clp.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form\cr
    \code{test <- initProbCLP()}.
}

\section{Slots}{
  \describe{
    \item{\code{clpPtrType}:}{
        Object of class \code{"character"} giving the pointer type.
    }
    \item{\code{clpPointer}:}{
        Object of class \code{"externalptr"} containig the pointer to a
        C structure.
    }
  }
}

\section{Methods}{
  \describe{
    \item{isCLPpointer}{\code{signature(object = "clpPtr")}:
        returns \code{TRUE} if \code{clpPointer(object)} is a pointer to a
        COIN-OR Clp problem object, otherwise \code{FALSE}.
    }
    \item{isNULLpointerCLP}{\code{signature(object = "clpPtr")}:
        returns \code{TRUE} if \code{clpPointer(object)} is a NULL pointer,
        otherwise \code{FALSE}.
    }
    \item{clpPointer}{\code{signature(object = "clpPtr")}:
        gets the \code{clpPointer} slot.
    }
    \item{clpPtrType}{\code{signature(object = "clpPtr")}:
        gets the \code{clpPtrType} slot.
    }
    \item{clpPtrType<-}{\code{signature(object = "clpPtr")}:
        sets the \code{clpPtrType} slot.
    }
  }
}

\references{
   The COIN-OR Clp home page at
   \url{https://projects.coin-or.org/Clp}
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Mayo Roettger <mayo.roettger@hhu.de>
}


\seealso{
  \code{\link{initProbCLP}}
}

\examples{
    showClass("clpPtr")
}

\keyword{ optimize }

