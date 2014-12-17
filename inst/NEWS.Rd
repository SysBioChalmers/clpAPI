%% clpAPI Version History

\name{NEWS}
\title{clpAPI News}
\encoding{UTF-8}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.6 2013-10-26}{
  \itemize{
    \item a minor bugfix in file \file{configure.ac}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.5 2013-09-08}{
  \itemize{
    \item removed argument \code{--enable-lapack} from file \file{configure.ac}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.4 2013-07-17}{
  \itemize{
    \item added function \code{probNameCLP()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.3 2013-05-23}{
  \itemize{
    \item prints the version number of COIN-OR Clp on startup
    \item added functions \code{dropNamesCLP()}, \code{copyNamesCLP()} and
          \code{lengthNamesCLP()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.2 2013-01-16}{
  \itemize{
    \item removed old file \file{NEWS.Rd} from top level directory
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.1 2012-04-21}{
  \itemize{
    \item bugfix in file \file{configure.ac} when using pkg-config
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.2.0 2012-04-20}{
  \itemize{
    \item updated \file{configure.ac} script

    \item file \file{NEWS} is now in \file{inst/NEWS.Rd}

    \item added finalizer for COIN-OR Clp problem objects

    \item Function \code{initProbCLP()} returns an instances of class
          \code{clpPtr}. This class has two data slots:
    \tabular{ll}{
      clpPtrType \tab a character string giving the pointer type and \cr
      clpPointer \tab a pointer to a C structure                     \cr
    }

    \item file \file{c2r.map} in \file{inst/}: mapping of filenames between
          C API and \pkg{clpAPI}, run
      \preformatted{
c2r <- system.file(package = "clpAPI", "c2r.map")
source(c2r)
      }
      in order to use original function names of the COIN-OR Clp C API.
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.7 2011-08-22}{
  \itemize{
    \item minor update in \file{configure.ac}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.6 2011-07-29}{
  \itemize{
    \item \file{configure.ac} script and \file{Makevars.win} are now compatible
          to COIN-OR Clp version 1.14.0
  }
}



% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.5 2011-07-28}{
  \itemize{
    \item added \file{NAMESPACE}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.4 2011-07-27}{
  \itemize{
    \item removed file \file{COPYING}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.3 2011-07-27}{
  \itemize{
    \item changes in \file{configure.ac}
    \item updated package vignette
    \item added files \file{INSTALL} and \file{README}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.2 2011-07-14}{
  \itemize{
    \item minor changes in \file{configure.ac}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.1 2011-07-12}{
  \itemize{
    \item renamed to \pkg{clpAPI}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.3 2011-05-12}{
  \itemize{
    \item added functions \code{getColDualCLP()}, \code{getRowPrimCLP()} and
          \code{getRowDualCLP()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.2 2010-12-15}{
  \itemize{
    \item fixed a bug in \file{configure.ac}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.1 2010-11-05}{
  \itemize{


    \item added functions \code{addRowsCLP()}, \code{addColsCLP()},
          \code{deleteRowsCLP()} and \code{deleteColsCLP()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.0 2010-08-19}{
  \itemize{
    \item initial public release
  }
}

