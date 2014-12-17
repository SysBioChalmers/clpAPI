/* init.c
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

#include <R.h>
#include <Rinternals.h>

#include "clpAPI.h"

#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"isCLPptr",                   (DL_FUNC) &isCLPptr,                   1},
    {"isNULLptr",                  (DL_FUNC) &isNULLptr,                  1},
    {"initCLP",                    (DL_FUNC) &initCLP,                    0},
    {"delProb",                    (DL_FUNC) &delProb,                    1},
    {"initProb",                   (DL_FUNC) &initProb,                   1},
    {"setObjDir",                  (DL_FUNC) &setObjDir,                  2},
    {"getObjDir",                  (DL_FUNC) &getObjDir,                  1},
    {"resize",                     (DL_FUNC) &resize,                     3},
    {"addRows",                    (DL_FUNC) &addRows,                    7},
    {"addCols",                    (DL_FUNC) &addCols,                    8},
    {"getNumRows",                 (DL_FUNC) &getNumRows,                 1},
    {"getNumCols",                 (DL_FUNC) &getNumCols,                 1},
    {"chgObjCoefs",                (DL_FUNC) &chgObjCoefs,                2},
    {"getObjCoefs",                (DL_FUNC) &getObjCoefs,                1},
    {"chgRowLower",                (DL_FUNC) &chgRowLower,                2},
    {"getRowLower",                (DL_FUNC) &getRowLower,                1},
    {"chgRowUpper",                (DL_FUNC) &chgRowUpper,                2},
    {"getRowUpper",                (DL_FUNC) &getRowUpper,                1},
    {"chgColLower",                (DL_FUNC) &chgColLower,                2},
    {"getColLower",                (DL_FUNC) &getColLower,                1},
    {"chgColUpper",                (DL_FUNC) &chgColUpper,                2},
    {"getColUpper",                (DL_FUNC) &getColUpper,                1},
    {"loadProblem",                (DL_FUNC) &loadProblem,               11},
    {"loadMatrix",                 (DL_FUNC) &loadMatrix,                 6},
    {"getNumNnz",                  (DL_FUNC) &getNumNnz,                  1},
    {"getVecStart",                (DL_FUNC) &getVecStart,                1},
    {"getInd",                     (DL_FUNC) &getInd,                     1},
    {"getVecLen",                  (DL_FUNC) &getVecLen,                  1},
    {"getNnz",                     (DL_FUNC) &getNnz,                     1},
    {"printModel",                 (DL_FUNC) &printModel,                 2},
    {"setLogLevel",                (DL_FUNC) &setLogLevel,                2},
    {"getLogLevel",                (DL_FUNC) &getLogLevel,                1},
    {"scaleModel",                 (DL_FUNC) &scaleModel,                 2},
    {"getScaleFlag",               (DL_FUNC) &getScaleFlag,               1},
    {"solveInitial",               (DL_FUNC) &solveInitial,               1},
    {"solveInitialDual",           (DL_FUNC) &solveInitialDual,           1},
    {"solveInitialPrimal",         (DL_FUNC) &solveInitialPrimal,         1},
    {"solveInitialBarrier",        (DL_FUNC) &solveInitialBarrier,        1},
    {"solveInitialBarrierNoCross", (DL_FUNC) &solveInitialBarrierNoCross, 1},
    {"dual",                       (DL_FUNC) &dual,                       2},
    {"primal",                     (DL_FUNC) &primal,                     2},
    {"idiot",                      (DL_FUNC) &idiot,                      2},
    {"getSolStatus",               (DL_FUNC) &getSolStatus,               1},
    {"getObjVal",                  (DL_FUNC) &getObjVal,                  1},
    {"getColPrim",                 (DL_FUNC) &getColPrim,                 1},
    {"getColDual",                 (DL_FUNC) &getColDual,                 1},
    {"getRowPrim",                 (DL_FUNC) &getRowPrim,                 1},
    {"getRowDual",                 (DL_FUNC) &getRowDual,                 1},
    {"delRows",                    (DL_FUNC) &delRows,                    3},
    {"delCols",                    (DL_FUNC) &delCols,                    3},
    {"readMPS",                    (DL_FUNC) &readMPS,                    4},
    {"saveModel",                  (DL_FUNC) &saveModel,                  2},
    {"restoreModel",               (DL_FUNC) &restoreModel,               2},
    {"version",                    (DL_FUNC) &version,                    0},
    {"dropNames",                  (DL_FUNC) &dropNames,                  1},
    {"copyNames",                  (DL_FUNC) &copyNames,                  3},
    {"lengthNames",                (DL_FUNC) &lengthNames,                1},
    {"rowName",                    (DL_FUNC) &rowName,                    3},
    {"colName",                    (DL_FUNC) &colName,                    3},
    {"probName",                   (DL_FUNC) &probName,                   3},
    {NULL, NULL, 0}
};


void R_init_clpAPI(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
