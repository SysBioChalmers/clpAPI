/* clpAPI.c
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


#include "clpAPI.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

static SEXP tagCLPprob;


/* -------------------------------------------------------------------------- */
/* Finalizer                                                                  */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* finalizer for clp problem objects */
static void clpProbFinalizer (SEXP lp) {
    if (!R_ExternalPtrAddr(lp)) {
        return;
    }
    else {
        delProb(lp);
    }
}


/* -------------------------------------------------------------------------- */
/* help functions                                                             */
/* -------------------------------------------------------------------------- */

/* check for pointer to clp */
SEXP isCLPptr(SEXP ptr) {

    SEXP out = R_NilValue;

    if ( (TYPEOF(ptr) == EXTPTRSXP) &&
         (R_ExternalPtrTag(ptr) == tagCLPprob) ) {
        out = Rf_ScalarLogical(1);
    }
    else {
        out = Rf_ScalarLogical(0);
    }

    return out;
}


/* check for NULL pointer */
SEXP isNULLptr(SEXP ptr) {

    SEXP out = R_NilValue;

    if ( (TYPEOF(ptr) == EXTPTRSXP) &&
         (R_ExternalPtrAddr(ptr) == NULL) ) {
        out = Rf_ScalarLogical(1);
    }
    else {
        out = Rf_ScalarLogical(0);
    }

    return out;
}


/* -------------------------------------------------------------------------- */
/* API-Functions                                                              */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* initialize clp */
SEXP initCLP(void) {
    tagCLPprob = Rf_install("TYPE_CLP_PROB");
    return R_NilValue;
}


/* -------------------------------------------------------------------------- */
/* remove problem object */
SEXP delProb(SEXP lp) {

    SEXP out = R_NilValue;
    Clp_Simplex *del = NULL;

    checkTypeOfProb(lp);

    del = R_ExternalPtrAddr(lp);

    Clp_deleteModel(del);
    R_ClearExternalPtr(lp);

    return out;
}


/* -------------------------------------------------------------------------- */
/* create new problem object */
SEXP initProb(SEXP ptrtype) {

    SEXP lpext = R_NilValue;
    SEXP ptr, class;

    Clp_Simplex *lp;

    /* create problem pointer */
    PROTECT(ptr = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(ptr, 0, STRING_ELT(ptrtype, 0));

    PROTECT(class = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, Rf_mkChar("clp_ptr"));

    lp = Clp_newModel();

    lpext = R_MakeExternalPtr(lp, tagCLPprob, R_NilValue);
    PROTECT(lpext);
    R_RegisterCFinalizerEx(lpext, clpProbFinalizer, TRUE);
    Rf_setAttrib(ptr, class, lpext);
    Rf_classgets(ptr, class);
    UNPROTECT(3);

    return ptr;
}


/* -------------------------------------------------------------------------- */
/* set optimization direction */
SEXP setObjDir(SEXP lp, SEXP dir) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_setOptimizationDirection(R_ExternalPtrAddr(lp), Rf_asReal(dir));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get optimization direction */
SEXP getObjDir(SEXP lp) {

    SEXP out = R_NilValue;
    double dir = 0;

    checkProb(lp);

    dir = Clp_optimizationDirection(R_ExternalPtrAddr(lp));

    out = Rf_ScalarReal(dir);

    return out;
}


/* -------------------------------------------------------------------------- */
/* resize the model */
SEXP resize(SEXP lp, SEXP nrows, SEXP ncols) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_resize(R_ExternalPtrAddr(lp), Rf_asInteger(nrows), Rf_asInteger(ncols));

    return out;
}


/* -------------------------------------------------------------------------- */
/* add rows */
SEXP addRows(SEXP lp, SEXP nrows,
             SEXP lb, SEXP ub, SEXP rowst, SEXP cols, SEXP val) {

    SEXP out = R_NilValue;

    const double *rlb  = REAL(lb);
    const double *rub  = REAL(ub);
    /* const int *rrowst  = INTEGER(rowst); */
    const CoinBigIndex *rrowst  = INTEGER(rowst);
    const int *rcols   = INTEGER(cols);
    const double *rval = REAL(val);

    checkProb(lp);

    Clp_addRows(R_ExternalPtrAddr(lp), Rf_asInteger(nrows),
                rlb, rub, rrowst, rcols, rval);

    return out;
}


/* -------------------------------------------------------------------------- */
/* add columns */
SEXP addCols(SEXP lp, SEXP ncols,
             SEXP lb, SEXP ub, SEXP obj, SEXP colst, SEXP rows, SEXP val) {

    SEXP out = R_NilValue;

    const double *rlb  = REAL(lb);
    const double *rub  = REAL(ub);
    const double *robj = REAL(obj);
    /* const int *rcolst  = INTEGER(colst); */
    const CoinBigIndex *rcolst  = INTEGER(colst);
    const int *rrows   = INTEGER(rows);
    const double *rval = REAL(val);

    checkProb(lp);

    Clp_addColumns(R_ExternalPtrAddr(lp), Rf_asInteger(ncols),
                   rlb, rub, robj, rcolst, rrows, rval);

    return out;
}
/* -------------------------------------------------------------------------- */
/* get maximum number of iterations */
SEXP getMaximumIterations(SEXP lp) {

    SEXP out = R_NilValue;
    int iterations = 0;

    checkProb(lp);

    iterations = maximumIterations(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(iterations);

    return out;
}

/* -------------------------------------------------------------------------- */
/* get maximum number of seconds */
SEXP getMaximumSeconds(SEXP lp) {

    SEXP out = R_NilValue;
    double seconds = 0;

    checkProb(lp);

    seconds = Clp_maximumSeconds(R_ExternalPtrAddr(lp));

    out = Rf_ScalarReal(seconds);

    return out;
}

/* -------------------------------------------------------------------------- */
/* get if maximum iteration (or time) bound was hit*/
SEXP getHitMaximumIterations(SEXP lp) {

    SEXP out = R_NilValue;

    checkProb(lp);

    if (Clp_hitMaximumIterations(R_ExternalPtrAddr(lp)))
        out = Rf_ScalarLogical(1);
    else
        out = Rf_ScalarLogical(0);
 
    return out;
}

/* -------------------------------------------------------------------------- */
/* get number of rows */
SEXP getNumRows(SEXP lp) {

    SEXP out = R_NilValue;
    int nrows = 0;

    checkProb(lp);

    nrows = Clp_numberRows(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(nrows);

    return out;
}


/* -------------------------------------------------------------------------- */
/* get number of columns */
SEXP getNumCols(SEXP lp) {

    SEXP out = R_NilValue;
    int ncols = 0;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ncols);

    return out;
}


/* -------------------------------------------------------------------------- */
/* set objective coefficients */
SEXP chgObjCoefs(SEXP lp, SEXP objCoef) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_chgObjCoefficients(R_ExternalPtrAddr(lp), REAL(objCoef));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get objective coefficients */
SEXP getObjCoefs(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    double *obj_coef;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));
    obj_coef = Clp_objective(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, ncols));
    for (k = 0; k < ncols; k++) {
        REAL(out)[k] = obj_coef[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* change row lower bounds */
SEXP chgRowLower(SEXP lp, SEXP rlb) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_chgRowLower(R_ExternalPtrAddr(lp), REAL(rlb));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get row lower bounds */
SEXP getRowLower(SEXP lp) {

    SEXP out = R_NilValue;

    int nrows, k;
    double *rlb;

    checkProb(lp);

    nrows = Clp_numberRows(R_ExternalPtrAddr(lp));
    rlb = Clp_rowLower(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, nrows));
    for (k = 0; k < nrows; k++) {
        REAL(out)[k] = rlb[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* change row upper bounds */
SEXP chgRowUpper(SEXP lp, SEXP rub) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_chgRowUpper(R_ExternalPtrAddr(lp), REAL(rub));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get row upper bounds */
SEXP getRowUpper(SEXP lp) {

    SEXP out = R_NilValue;

    int nrows, k;
    double *rub;

    checkProb(lp);

    nrows = Clp_numberRows(R_ExternalPtrAddr(lp));
    rub = Clp_rowUpper(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, nrows));
    for (k = 0; k < nrows; k++) {
        REAL(out)[k] = rub[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* change column lower bounds */
SEXP chgColLower(SEXP lp, SEXP lb) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_chgColumnLower(R_ExternalPtrAddr(lp), REAL(lb));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get column lower bounds */
SEXP getColLower(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    double *lb;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));
    lb = Clp_columnLower(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, ncols));
    for (k = 0; k < ncols; k++) {
        REAL(out)[k] = lb[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* change column upper bounds */
SEXP chgColUpper(SEXP lp, SEXP ub) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_chgColumnUpper(R_ExternalPtrAddr(lp), REAL(ub));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get column upper bounds */
SEXP getColUpper(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    double *ub;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));
    ub = Clp_columnUpper(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, ncols));
    for (k = 0; k < ncols; k++) {
        REAL(out)[k] = ub[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* load a complete Problem or at least the constraint matrix */
SEXP loadProblem(SEXP lp, SEXP ncols, SEXP nrows, SEXP ia, SEXP ja, SEXP ra,
                 SEXP clb, SEXP cub, SEXP obj_coef, SEXP rlb, SEXP rub) {

    SEXP out = R_NilValue;

    int *ria = INTEGER(ia);
    /* int *rja = INTEGER(ja); */
    CoinBigIndex *rja = INTEGER(ja);
    double *rra = REAL(ra);
    double *rclb;
    double *rcub;
    double *robj_coef;
    double *rrlb;
    double *rrub;
/*     double *rclb = REAL(clb); */
/*     double *rcub = REAL(cub); */
/*     double *robj_coef = REAL(obj_coef); */
/*     double *rrlb = REAL(rlb); */
/*     double *rrub = REAL(rub); */

    checkProb(lp);

    if (clb == R_NilValue) {
        rclb = NULL;
    }
    else {
        rclb = REAL(clb);
    }

    if (cub == R_NilValue) {
        rcub = NULL;
    }
    else {
        rcub = REAL(cub);
    }

    if (obj_coef == R_NilValue) {
        robj_coef = NULL;
    }
    else {
        robj_coef = REAL(obj_coef);
    }

    if (rlb == R_NilValue) {
        rrlb = NULL;
    }
    else {
        rrlb = REAL(rlb);
    }

    if (rub == R_NilValue) {
        rrub = NULL;
    }
    else {
        rrub = REAL(rub);
    }

    Clp_loadProblem(R_ExternalPtrAddr(lp), Rf_asInteger(ncols),
                    Rf_asInteger(nrows), rja, ria, rra, rclb, rcub,
                    robj_coef, rrlb, rrub
                   );

    return out;
}


/* -------------------------------------------------------------------------- */
/* load a complete Problem or at least the constraint matrix */
SEXP loadMatrix(SEXP lp, SEXP ncols, SEXP nrows, SEXP ia, SEXP ja, SEXP ra) {

    SEXP out = R_NilValue;

    int *ria = INTEGER(ia);
    /* int *rja = INTEGER(ja); */
    CoinBigIndex *rja = INTEGER(ja);
    double *rra = REAL(ra);

    checkProb(lp);

    Clp_loadProblem(R_ExternalPtrAddr(lp),
                    Rf_asInteger(ncols), Rf_asInteger(nrows),
                    rja, ria, rra, NULL, NULL, NULL, NULL, NULL
                   );

    return out;
}


/* -------------------------------------------------------------------------- */
/* get number of non zero elements in the contraint matrix */
SEXP getNumNnz(SEXP lp) {

    SEXP out = R_NilValue;
    CoinBigIndex nnz;

    checkProb(lp);

    nnz = Clp_getNumElements(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(nnz);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Column starts in constraint matrix (ja(-a) in column major order format) */
SEXP getVecStart(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    /* const int *vec_start; */
    const CoinBigIndex *vec_start;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp)) + 1;
    vec_start = Clp_getVectorStarts(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(INTSXP, ncols));
    for (k = 0; k < ncols; k++) {
        INTEGER(out)[k] = vec_start[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Row indices in the constraint matrix (ia(-1) in column major order format) */
SEXP getInd(SEXP lp) {

    SEXP out = R_NilValue;

    int nnz, k;
    const int *index;

    checkProb(lp);

    nnz = Clp_getNumElements(R_ExternalPtrAddr(lp));
    index = Clp_getIndices(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(INTSXP, nnz));
    for (k = 0; k < nnz; k++) {
        INTEGER(out)[k] = index[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Vector (number of nnz per column) length in the constraint matrix
   (lg in column major order format) */
SEXP getVecLen(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    const int *vec_len;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));
    vec_len = Clp_getVectorLengths(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(INTSXP, ncols));
    for (k = 0; k < ncols; k++) {
        INTEGER(out)[k] = vec_len[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Non zero elements in the constraint matrix (ar in column major format) */
SEXP getNnz(SEXP lp) {

    SEXP out = R_NilValue;

    int nnz, k;
    const double *n_elem;

    checkProb(lp);

    nnz = Clp_getNumElements(R_ExternalPtrAddr(lp));
    n_elem = Clp_getElements(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, nnz));
    for (k = 0; k < nnz; k++) {
        REAL(out)[k] = n_elem[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* print model */
SEXP printModel(SEXP lp, SEXP prefix) {

    SEXP out = R_NilValue;
    const char *rprefix = CHAR(STRING_ELT(prefix, 0));

    checkProb(lp);

    Clp_printModel(R_ExternalPtrAddr(lp), rprefix);
    /* Clp_printModel(R_ExternalPtrAddr(lp), "CLPmodel"); */

    return out;
}

/* -------------------------------------------------------------------------- */
/* set number of iterations */
SEXP setNumberIterations(SEXP lp, SEXP iterations) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_setNumberIterations(R_ExternalPtrAddr(lp), Rf_asInteger(iterations));

    return out;
}
/* -------------------------------------------------------------------------- */
/* set maximal number of iterations */
SEXP setMaximumIterations(SEXP lp, SEXP iterations) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_setMaximumIterations(R_ExternalPtrAddr(lp), Rf_asInteger(iterations));

    return out;
}
/* -------------------------------------------------------------------------- */
/* set maximal duration in seconds */
SEXP setMaximumSeconds(SEXP lp, SEXP seconds) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_setMaximumSeconds(R_ExternalPtrAddr(lp), Rf_asReal(seconds));

    return out;
}

/* -------------------------------------------------------------------------- */
/* amount of print out */
SEXP setLogLevel(SEXP lp, SEXP amount) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_setLogLevel(R_ExternalPtrAddr(lp), Rf_asInteger(amount));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get log level */
SEXP getLogLevel(SEXP lp) {

    SEXP out = R_NilValue;
    int amount;

    checkProb(lp);

    amount = Clp_logLevel(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(amount);

    return out;
}


/* -------------------------------------------------------------------------- */
/* set or unset scaling */
SEXP scaleModel(SEXP lp, SEXP mode) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_scaling(R_ExternalPtrAddr(lp), Rf_asInteger(mode));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get scale flag */
SEXP getScaleFlag(SEXP lp) {

    SEXP out = R_NilValue;
    int flag;

    checkProb(lp);

    flag = Clp_scalingFlag(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(flag);

    return out;
}


/* -------------------------------------------------------------------------- */
/* solve model with general solve algorithm */
SEXP solveInitial(SEXP lp) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_initialSolve(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Dual initial solve */
SEXP solveInitialDual(SEXP lp) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_initialDualSolve(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Primal initial solve */
SEXP solveInitialPrimal(SEXP lp) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_initialPrimalSolve(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Barrier initial solve */
SEXP solveInitialBarrier(SEXP lp) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_initialBarrierSolve(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Barrier initial solve, no crossover */
SEXP solveInitialBarrierNoCross(SEXP lp) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_initialBarrierNoCrossSolve(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Dual slgorithm */
SEXP dual(SEXP lp, SEXP ifValP) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_dual(R_ExternalPtrAddr(lp), Rf_asInteger(ifValP));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* Dual slgorithm */
SEXP primal(SEXP lp, SEXP ifValP) {

    SEXP out = R_NilValue;
    int ret;

    checkProb(lp);

    ret = Clp_primal(R_ExternalPtrAddr(lp), Rf_asInteger(ifValP));

    out = Rf_ScalarInteger(ret);

    return out;
}


/* -------------------------------------------------------------------------- */
/* solve problem using the idiot code */
SEXP idiot(SEXP lp, SEXP thd) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_idiot(R_ExternalPtrAddr(lp), Rf_asInteger(thd));

    return out;
}


/* -------------------------------------------------------------------------- */
/* get solution status */
SEXP getSolStatus(SEXP lp) {

    SEXP out = R_NilValue;
    int stat;

    checkProb(lp);

    stat = Clp_status(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(stat);

    return out;
}


/* -------------------------------------------------------------------------- */
/* solve model with general solve algorithm */
SEXP getObjVal(SEXP lp) {

    SEXP out = R_NilValue;
    double obj;

    checkProb(lp);

    obj = Clp_objectiveValue(R_ExternalPtrAddr(lp));

    out = Rf_ScalarReal(obj);

    return out;
}


/* -------------------------------------------------------------------------- */
/* get column primal solution */
SEXP getColPrim(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    double *col_prim;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));
    col_prim = Clp_primalColumnSolution(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, ncols));
    for (k = 0; k < ncols; k++) {
        REAL(out)[k] = col_prim[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* get column dual solution */
SEXP getColDual(SEXP lp) {

    SEXP out = R_NilValue;

    int ncols, k;
    double *col_dual;

    checkProb(lp);

    ncols = Clp_numberColumns(R_ExternalPtrAddr(lp));
    col_dual = Clp_dualColumnSolution(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, ncols));
    for (k = 0; k < ncols; k++) {
        REAL(out)[k] = col_dual[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* get row primal solution */
SEXP getRowPrim(SEXP lp) {

    SEXP out = R_NilValue;

    int nrows, k;
    double *row_prim;

    checkProb(lp);

    nrows = Clp_numberRows(R_ExternalPtrAddr(lp));
    row_prim = Clp_primalRowSolution(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, nrows));
    for (k = 0; k < nrows; k++) {
        REAL(out)[k] = row_prim[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* get row dual solution */
SEXP getRowDual(SEXP lp) {

    SEXP out = R_NilValue;

    int nrows, k;
    double *row_dual;

    checkProb(lp);

    nrows = Clp_numberRows(R_ExternalPtrAddr(lp));
    row_dual = Clp_dualRowSolution(R_ExternalPtrAddr(lp));

    PROTECT(out = Rf_allocVector(REALSXP, nrows));
    for (k = 0; k < nrows; k++) {
        REAL(out)[k] = row_dual[k];
    }
    UNPROTECT(1);

    return out;
}


/* -------------------------------------------------------------------------- */
/* delete rows */
SEXP delRows(SEXP lp, SEXP num, SEXP i) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_deleteRows(R_ExternalPtrAddr(lp), Rf_asInteger(num), INTEGER(i));

    return out;
}


/* -------------------------------------------------------------------------- */
/* delete columns */
SEXP delCols(SEXP lp, SEXP num, SEXP j) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_deleteColumns(R_ExternalPtrAddr(lp), Rf_asInteger(num), INTEGER(j));

    return out;
}


/* -------------------------------------------------------------------------- */
/* read problem data in MPS format */
SEXP readMPS(SEXP lp, SEXP fname, SEXP keepNames, SEXP ignoreErrors) {

    SEXP out = R_NilValue;
    const char *rfname = CHAR(STRING_ELT(fname, 0));
    int check = 0;

    checkProb(lp);

    check = Clp_readMps(R_ExternalPtrAddr(lp), rfname, Rf_asInteger(keepNames),
                        Rf_asInteger(ignoreErrors)
                       );

    out = Rf_ScalarInteger(check);

    return out;
}


/* -------------------------------------------------------------------------- */
/* save model to file */
SEXP saveModel(SEXP lp, SEXP fname) {

    SEXP out = R_NilValue;
    const char *rfname = CHAR(STRING_ELT(fname, 0));
    int check = 0;

    checkProb(lp);

    check = Clp_saveModel(R_ExternalPtrAddr(lp), rfname);

    out = Rf_ScalarInteger(check);

    return out;
}


/* -------------------------------------------------------------------------- */
/* restore model from file */
SEXP restoreModel(SEXP lp, SEXP fname) {

    SEXP out = R_NilValue;
    const char *rfname = CHAR(STRING_ELT(fname, 0));
    int check = 0;

    checkProb(lp);

    check = Clp_restoreModel(R_ExternalPtrAddr(lp), rfname);

    out = Rf_ScalarInteger(check);

    return out;
}


/* -------------------------------------------------------------------------- */
/* get COIN OR Clp version */
SEXP version() {

    SEXP out = R_NilValue;
    
    out = Rf_mkString(CLP_VERSION);

    return out;
}


/* -------------------------------------------------------------------------- */
/* drop names */
SEXP dropNames(SEXP lp) {

    SEXP out = R_NilValue;

    checkProb(lp);

    Clp_dropNames(R_ExternalPtrAddr(lp));

    return out;
}


/* -------------------------------------------------------------------------- */
/* copy in names */
SEXP copyNames(SEXP lp, SEXP cnames, SEXP rnames) {

    SEXP out = R_NilValue;

    int k, numcn, numrn;
    const char **rcnames;
    const char ** rrnames;

    checkProb(lp);

    numcn = Rf_length(cnames);
    rcnames = R_Calloc(numcn, const char *);
    for (k = 0; k < numcn; k++) {
        rcnames[k] = CHAR(STRING_ELT(cnames, k));
    }

    numrn = Rf_length(rnames);
    rrnames = R_Calloc(numrn, const char *);
    for (k = 0; k < numrn; k++) {
        rrnames[k] = CHAR(STRING_ELT(rnames, k));
    }

    Clp_copyNames(R_ExternalPtrAddr(lp), rrnames, rcnames);

    if (cnames != R_NilValue) {
        R_Free(rcnames);
    }
    if (rnames != R_NilValue) {
        R_Free(rrnames);
    }

    return out;
}


/* -------------------------------------------------------------------------- */
/* length of names */
SEXP lengthNames(SEXP lp) {

    SEXP out = R_NilValue;
    int ncnames;

    checkProb(lp);

    ncnames = Clp_lengthNames(R_ExternalPtrAddr(lp));

    out = Rf_ScalarInteger(ncnames);

    return out;
}


/* -------------------------------------------------------------------------- */
/* fill in row name */
SEXP rowName(SEXP lp, SEXP i, SEXP rname) {

    SEXP out = R_NilValue;

    const char *rrname = CHAR(STRING_ELT(rname, 0));

    checkProb(lp);

    Clp_rowName(R_ExternalPtrAddr(lp), Rf_asInteger(i), (char *) rrname);

    return out;
}


/* -------------------------------------------------------------------------- */
/* fill in column name */
SEXP colName(SEXP lp, SEXP j, SEXP cname) {

    SEXP out = R_NilValue;

    const char *rcname = CHAR(STRING_ELT(cname, 0));

    checkProb(lp);

    Clp_columnName(R_ExternalPtrAddr(lp), Rf_asInteger(j), (char *) rcname);

    return out;
}


/* -------------------------------------------------------------------------- */
/* fill in problem name */
SEXP probName(SEXP lp, SEXP nc, SEXP pname) {

    SEXP out = R_NilValue;

    const char *rpname = CHAR(STRING_ELT(pname, 0));

    checkProb(lp);

    Clp_problemName(R_ExternalPtrAddr(lp), Rf_asInteger(nc), (char *) rpname);

    return out;
}



#ifdef HAVE_CLP_EXT1_17_2

/* NEW in Clp-1.17.2 */
/* -------------------------------------------------------------------------- */
/* fill in row name */
SEXP setRowName(SEXP lp, SEXP i, SEXP rname) {

    SEXP out = R_NilValue;

    const char *rrname = CHAR(STRING_ELT(rname, 0));

    checkProb(lp);

    Clp_setRowName(R_ExternalPtrAddr(lp), Rf_asInteger(i), (char *) rrname);

    return out;
}


/* NEW in Clp-1.17.2 */
/* -------------------------------------------------------------------------- */
/* fill in column name */
SEXP setColName(SEXP lp, SEXP j, SEXP cname) {

    SEXP out = R_NilValue;

    const char *rcname = CHAR(STRING_ELT(cname, 0));

    checkProb(lp);

    Clp_setColumnName(R_ExternalPtrAddr(lp), Rf_asInteger(j), (char *) rcname);

    return out;
}


/* NEW in Clp-1.17.2 */
/* -------------------------------------------------------------------------- */
/* Write an mps file to the given filename */
SEXP writeMps(SEXP lp, SEXP filename, SEXP formatType, SEXP numberAcross, SEXP objSense) {
    
    int check = 0;
    
    const char *rfilename = CHAR(STRING_ELT(filename, 0));
    
    checkProb(lp);
    
    check = Clp_writeMps(R_ExternalPtrAddr(lp), rfilename, Rf_asInteger(formatType), Rf_asInteger(numberAcross), Rf_asReal(objSense));
    
    return Rf_ScalarInteger(check);
}


/* NEW in Clp-1.17.2 */
/* -------------------------------------------------------------------------- */
/* Change matrix coefficients */
SEXP modifyCoefficient(SEXP lp, SEXP row, SEXP column, SEXP newElement, SEXP keepZero) {
    
    SEXP out = R_NilValue;
    
    bool rkeepZero = Rf_asLogical(keepZero);
    
    checkProb(lp);
    
    Clp_modifyCoefficient(R_ExternalPtrAddr(lp), Rf_asInteger(row), Rf_asInteger(column), Rf_asReal(newElement), rkeepZero);
    
    return out;
}

/* -------------------------------------------------------------------------- */
/* check for functionality of new clp functions */
SEXP isAvailableFunc(SEXP funcname) {
    
    SEXP out = R_NilValue;
    
    const char *rfuncname = CHAR(STRING_ELT(funcname, 0));
    
    if (strcmp(rfuncname,"setRowNameCLP") == 0) {
        out = Rf_ScalarLogical(1);
    } else if (strcmp(rfuncname,"setColNameCLP") == 0) {
        out = Rf_ScalarLogical(1);
    } else if (strcmp(rfuncname,"writeMpsCLP") == 0) {
        out = Rf_ScalarLogical(1);
    } else if (strcmp(rfuncname,"modifyCoefficientCLP") == 0) {
        out = Rf_ScalarLogical(1);
    }
    
    return out;
}

#else /* not CLP_EXT1_17_2 */

/* dummy function */
SEXP setRowName(SEXP lp, SEXP i, SEXP rname) {
    SEXP out = R_NilValue;
    return out;
}

/* dummy function */
SEXP setColName(SEXP lp, SEXP j, SEXP cname) {
    SEXP out = R_NilValue;
    return out;
}

/* dummy function */
SEXP writeMps(SEXP lp, SEXP filename, SEXP formatType, SEXP numberAcross, SEXP objSense) {
    SEXP out = Rf_ScalarInteger(1);
    return out;
}

/* dummy function */
SEXP modifyCoefficient(SEXP lp, SEXP row, SEXP column, SEXP newElement, SEXP keepZero) {
    SEXP out = R_NilValue;
    return out;
}

/* check for functionality of new clp functions */
SEXP isAvailableFunc(SEXP funcname) {
    
    SEXP out = R_NilValue;
    
    const char *rfuncname = CHAR(STRING_ELT(funcname, 0));
    
    if (strcmp(rfuncname,"setRowNameCLP") == 0) {
        out = Rf_ScalarLogical(0);
    } else if (strcmp(rfuncname,"setColNameCLP") == 0) {
        out = Rf_ScalarLogical(0);
    } else if (strcmp(rfuncname,"writeMpsCLP") == 0) {
        out = Rf_ScalarLogical(0);
    } else if (strcmp(rfuncname,"modifyCoefficientCLP") == 0) {
        out = Rf_ScalarLogical(0);
    }
    
    return out;
}

#endif /* HAVE_CLP_EXT1_17_2 */
