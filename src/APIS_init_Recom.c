#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>


/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void get_individual_mendelian_probability_2n(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_individual_mendelian_probability_3n(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void*);
extern void get_mendelian_probability_2n(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void get_mendelian_probability_3n(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"get_individual_mendelian_probability_2n", (DL_FUNC) &get_individual_mendelian_probability_2n, 10},
    {"get_individual_mendelian_probability_3n", (DL_FUNC) &get_individual_mendelian_probability_3n, 11},
    {"get_mendelian_probability_2n",            (DL_FUNC) &get_mendelian_probability_2n,            11},
    {"get_mendelian_probability_3n",            (DL_FUNC) &get_mendelian_probability_3n,            11},
    {NULL, NULL, 0}
};


void R_init_APIS(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
