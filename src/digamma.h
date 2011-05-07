#define MATHLIB_STANDALONE TRUE
#include <Rmath.h>

#define DIGAMMA_TABLE_SIZE 1000000
#define DIGAMMA_TABLE_MIN 0 // don't change this! even though 1 is min necessary.
#define DIGAMMA_TABLE_MID 10
#define DIGAMMA_TABLE_MAX 50
#define ROUND(x) round(x)

void *make_digamma_table() ;
double digamma_lookup(double x) ;
void digamma_test(double *x, int *n, double *ans) ;
FILE *fp ;

