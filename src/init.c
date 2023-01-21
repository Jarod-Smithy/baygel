#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _baygel_blockBAGR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBSGR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_mvrnormArma(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_baygel_blockBAGR",   (DL_FUNC) &_baygel_blockBAGR,   6},
  {"_baygel_blockBSGR",   (DL_FUNC) &_baygel_blockBSGR,   6},
  {"_baygel_mvrnormArma", (DL_FUNC) &_baygel_mvrnormArma, 3},
  {NULL, NULL, 0}
};

void R_init_baygel(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
