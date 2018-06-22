/* clpAPI.h
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


#include "clpR.h"


/* -------------------------------------------------------------------------- */
/* help functions                                                             */
/* -------------------------------------------------------------------------- */

/* check for pointer to clp */
SEXP isCLPptr(SEXP ptr);

/* check for NULL pointer */
SEXP isNULLptr(SEXP ptr);


/* -------------------------------------------------------------------------- */
/* API functions                                                              */
/* -------------------------------------------------------------------------- */

/* initialize clp */
SEXP initCLP(void);

/* remove problem object */
SEXP delProb(SEXP lp);

/* create new problem object */
SEXP initProb(SEXP ptrtype);

/* set optimization direction */
SEXP setObjDir(SEXP lp, SEXP dir);

/* get optimization direction */
SEXP getObjDir(SEXP lp);

/* resize the model */
SEXP resize(SEXP lp, SEXP nrows, SEXP ncols);

/* add rows */
SEXP addRows(SEXP lp, SEXP nrows,
             SEXP lb, SEXP ub, SEXP rowst, SEXP cols, SEXP val);

/* add columns */
SEXP addCols(SEXP lp, SEXP ncols,
             SEXP lb, SEXP ub, SEXP obj, SEXP colst, SEXP rows, SEXP val);

/* get maximum number of iterations */
SEXP getMaximumIterations(SEXP lp);

/* get maximum number of seconds */
SEXP getMaximumSeconds(SEXP lp);

/* get if maxium iteration bound was hit*/
SEXP getHitMaximumIterations(SEXP lp);

/* get number of rows */
SEXP getNumRows(SEXP lp);

/* get number of columns */
SEXP getNumCols(SEXP lp);

/* set objective coefficients */
SEXP chgObjCoefs(SEXP lp, SEXP objCoef);

/* get objective coefficients */
SEXP getObjCoefs(SEXP lp);

/* change row lower bounds */
SEXP chgRowLower(SEXP lp, SEXP rlb);

/* get row lower bounds */
SEXP getRowLower(SEXP lp);

/* change row upper bounds */
SEXP chgRowUpper(SEXP lp, SEXP rub);

/* get row upper bounds */
SEXP getRowUpper(SEXP lp);

/* change column lower bounds */
SEXP chgColLower(SEXP lp, SEXP lb);

/* get column lower bounds */
SEXP getColLower(SEXP lp);

/* change column upper bounds */
SEXP chgColUpper(SEXP lp, SEXP ub);

/* get column upper bounds */
SEXP getColUpper(SEXP lp);

/* load a complete Problem or at least the constraint matrix */
SEXP loadProblem(SEXP lp, SEXP ncols, SEXP nrows, SEXP ia, SEXP ja, SEXP ra,
                 SEXP clb, SEXP cub, SEXP obj_coef, SEXP rlb, SEXP rub);

/* load a complete Problem or at least the constraint matrix */
SEXP loadMatrix(SEXP lp, SEXP ncols, SEXP nrows, SEXP ia, SEXP ja, SEXP ra);

/* get number of non zero elements in the contraint matrix */
SEXP getNumNnz(SEXP lp);

/* Column starts in constraint matrix (ja(-a) in column major order format) */
SEXP getVecStart(SEXP lp);
 
/* Row indices in the constraint matrix (ia(-1) in column major order format) */
SEXP getInd(SEXP lp);
 
/* Vector (number of nnz per column) length in the constraint matrix
   (lg in column major order format) */
SEXP getVecLen(SEXP lp);
 
/* Non zero elements in the constraint matrix (ar in column major format) */
SEXP getNnz(SEXP lp);

/* print model */
SEXP printModel(SEXP lp, SEXP prefix);

/* set number of iterations */
SEXP setNumberIterations(SEXP lp, SEXP iterations);

/* set maximal number of iterations */
SEXP setMaximumIterations(SEXP lp, SEXP iterations);

/* set maximal duration in seconds */
SEXP setMaximumSeconds(SEXP lp, SEXP seconds);

/* amount of print out */
SEXP setLogLevel(SEXP lp, SEXP amount);

/* get log level */
SEXP getLogLevel(SEXP lp);

/* set or unset scaling */
SEXP scaleModel(SEXP lp, SEXP mode);

/* get scale flag */
SEXP getScaleFlag(SEXP lp);

/* solve model with general solve algorithm */
SEXP solveInitial(SEXP lp);

/* Dual initial solve */
SEXP solveInitialDual(SEXP lp);

/* Primal initial solve */
SEXP solveInitialPrimal(SEXP lp);

/* Barrier initial solve */
SEXP solveInitialBarrier(SEXP lp);

/* Barrier initial solve, no crossover */
SEXP solveInitialBarrierNoCross(SEXP lp);

/* Dual slgorithm */
SEXP dual(SEXP lp, SEXP ifValP);

/* Dual slgorithm */
SEXP primal(SEXP lp, SEXP ifValP);

/* solve problem using the idiot code */
SEXP idiot(SEXP lp, SEXP thd);

/* get solution status */
SEXP getSolStatus(SEXP lp);

/* solve model with general solve algorithm */
SEXP getObjVal(SEXP lp);

/* get column primal solution */
SEXP getColPrim(SEXP lp);

/* get column dual solution */
SEXP getColDual(SEXP lp);

/* get row primal solution */
SEXP getRowPrim(SEXP lp);

/* get row dual solution */
SEXP getRowDual(SEXP lp);

/* delete rows */
SEXP delRows(SEXP lp, SEXP num, SEXP i);

/* delete columns */
SEXP delCols(SEXP lp, SEXP num, SEXP j);

/* read problem data in MPS format */
SEXP readMPS(SEXP lp, SEXP fname, SEXP keepNames, SEXP ignoreErrors);

/* save model to file */
SEXP saveModel(SEXP lp, SEXP fname);

/* restore model from file */
SEXP restoreModel(SEXP lp, SEXP fname);

/* get COIN OR Clp version */
SEXP version();

/* drop names */
SEXP dropNames(SEXP lp);

/* copy in names */
SEXP copyNames(SEXP lp, SEXP cnames, SEXP rnames);

/* length of names */
SEXP lengthNames(SEXP lp);

/* fill in row name */
SEXP rowName(SEXP lp, SEXP i, SEXP rname);

/* fill in column name */
SEXP colName(SEXP lp, SEXP j, SEXP cname);

/* fill in problem name */
SEXP probName(SEXP lp, SEXP nc, SEXP pname);

