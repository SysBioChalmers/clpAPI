/* clpR.h
   R interface to COIN-OR Clp.
 
   Copyright (C) 2011-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
   Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
   All right reserved.
   Email: geliudie@uni-duesseldorf.de
 
   This file is part of clpAPI.
 
   ClpAPI is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   ClpAPI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with clpAPI.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>

/* include config.h definitions */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

/* include stdbool.h */
#ifdef HAVE_STDBOOL_H 
#include <stdbool.h>
#endif /* HAVE_STDBOOL_H */

#include <Clp_C_Interface.h>
#include <ClpConfig.h>

/* avoid remapping of Rf_<function> to <function> in R header files */
#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif /* R_NO_REMAP */

/* use strict R headers */
#ifndef STRICT_R_HEADERS
#define STRICT_R_HEADERS
#endif /* STRICT_R_HEADERS */

#include <R.h>
#include <Rinternals.h>


/* -------------------------------------------------------------------------- */
/* NULL */
#define checkIfNil(cp) do { \
    if (R_ExternalPtrAddr(cp) == NULL) \
        Rf_error("You passed a nil value!"); \
} while (0)


/* -------------------------------------------------------------------------- */
/* problem */
#define checkTypeOfProb(cp) do { \
    if ( (TYPEOF(cp) != EXTPTRSXP) || (R_ExternalPtrTag(cp) != tagCLPprob) ) \
        Rf_error("You must pass a clp problem!"); \
} while (0)

#define checkProb(p) do { \
    checkIfNil(p); \
    checkTypeOfProb(p); \
} while (0)

