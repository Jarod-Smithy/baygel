#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _baygel_blockBAGENI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBAGENII(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBAGL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBAGR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBGEN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBGL(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _baygel_blockBGR(SEXP, SEXP, SEXP, SEXP, SEXP);



extern SEXP _baygel_mvrnormArma(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_baygel_blockBAGENI",  (DL_FUNC) &_baygel_blockBAGENI,  8},
  {"_baygel_blockBAGENII", (DL_FUNC) &_baygel_blockBAGENII, 6},
  {"_baygel_blockBAGL",    (DL_FUNC) &_baygel_blockBAGL,    6},
  {"_baygel_blockBAGR",    (DL_FUNC) &_baygel_blockBAGR,    6},
  {"_baygel_blockBGEN",    (DL_FUNC) &_baygel_blockBGEN,    6},
  {"_baygel_blockBGL",     (DL_FUNC) &_baygel_blockBGL,     5},
  {"_baygel_blockBGR",     (DL_FUNC) &_baygel_blockBGR,     5},
  {NULL, NULL, 0}
};

void R_init_baygel(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
