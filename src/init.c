#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _baygel_ABGR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBSGR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_mvrnormArma(SEXP, SEXP, SEXP);

