#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Fortran calls */
  extern void F77_NAME(cormarkovchain)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(disag3daygwexprec_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"cormarkovchain",      (DL_FUNC) &F77_NAME(cormarkovchain),       6},
  {"disag3daygwexprec_f", (DL_FUNC) &F77_NAME(disag3daygwexprec_f), 13},
  {NULL, NULL, 0}
};

void R_init_GWEX(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}